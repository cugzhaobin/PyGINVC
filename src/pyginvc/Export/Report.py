#!/usr/bin/env python
"""HtmlReport -- build a self-contained HTML report of an inversion run.

The report collects run metadata, key metrics, result tables, figures and
download links, then renders them through a single HTML template.  It is
independent of Fault/GeoData/solution classes: ``Export.Output`` fills it in,
but any caller with paths and numbers can use it too.

Typical use::

    rep = HtmlReport(outdir, title='Coseismic slip inversion')
    rep.add_meta('Run', '2026-11-14 10:00')
    rep.add_metric('GPS stations', 128)
    rep.add_table('Inversion summary', ['smoothing', 'WRSS'], rows)
    rep.add_figures_from_dir()
    rep.write()

Set ``embed=True`` to inline every image as base64 so the HTML can be shared
as a single file.
"""
import base64
import glob
import html
import logging
import os

try:
    import numpy as np
except Exception:  # pragma: no cover - numpy is a hard dependency anyway
    np = None

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S')

_TEMPLATE_CACHE = None

_IMAGE_MIME = {
    '.png': 'image/png',
    '.jpg': 'image/jpeg',
    '.jpeg': 'image/jpeg',
    '.gif': 'image/gif',
    '.svg': 'image/svg+xml',
    '.webp': 'image/webp',
}

# Figures produced by Export.View, in the order they should appear.
FIGURE_PATTERNS = (
    ('GPS observation vs prediction', 'gps_mod_obs*.png'),
    ('InSAR observation / prediction / residual', 'sar_mod_obs*.png'),
    ('InSAR ramp', 'sar_ramp*.png'),
    ('3D slip distribution', 'slip_3d*.png'),
    ('Trade-off (L-) curve', 'lcurve*.png'),
    ('MCMC traces', 'mcmc_*trace*.png'),
    ('MCMC corner plot', 'mcmc_corner*.png'),
)

# Result files offered as downloads, in the order they should appear.
DOWNLOAD_PATTERNS = (
    'solutions.h5', 'mcmc_solutions.h5', 'mcmc_chain.h5',
    'summary', 'ruff_misfit',
    'slipmodel*.faultgeom', 'slipmodel*.tri',
    'ss_slip_llh.gmtlin', 'ds_slip_llh.gmtlin', 'ts_slip_llh.gmtlin',
    'gps*.gmtvec', 'sar*.xyz', 'lev*.xyz',
)


def _load_template():
    '''Read the HTML template once, falling back to an inline copy.'''
    global _TEMPLATE_CACHE
    if _TEMPLATE_CACHE is not None:
        return _TEMPLATE_CACHE
    try:
        from importlib.resources import files
        _TEMPLATE_CACHE = files('pyginvc.Export').joinpath(
            'templates/report.html').read_text(encoding='utf-8')
    except Exception:
        _TEMPLATE_CACHE = (
            '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">'
            '<title>$title</title><style>'
            'body{font-family:sans-serif;margin:24px;background:#f6f8fa}'
            'figure{margin:0 0 24px;background:#fff;border:1px solid #ddd;'
            'border-radius:8px;padding:10px}img{max-width:100%;height:auto}'
            'table{border-collapse:collapse;font-size:13px;background:#fff}'
            'th,td{padding:6px 10px;border-bottom:1px solid #eee}'
            '</style></head><body><h1>$title</h1><div>$subtitle</div>$meta'
            '$body</body></html>')
    return _TEMPLATE_CACHE


def _fmt(value):
    '''Format a scalar for display inside a table or card.'''
    if np is not None and isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if value == 0:
            return '0'
        a = abs(value)
        if a >= 1e5 or a < 1e-3:
            return '%.3e' % value
        return '%.4g' % value
    return html.escape(str(value))


class HtmlReport(object):
    '''
    Collect pieces of a run and write ``report.html``.

    Input:
        outdir : string
            Directory that holds the figures/result files; the report is
            written here and links are relative to it.
        title  : string
            Page title.
        embed  : bool
            Inline images as base64 when True (single-file report).
    '''

    def __init__(self, outdir='.', title='PyGINVC report', embed=False):
        self.outdir   = os.path.abspath(outdir or '.')
        self.title    = title
        self.embed    = embed
        self._meta    = []       # [(label, value)]
        self._metrics = []       # [(label, value, unit)]
        self._items   = []       # ordered blocks
        self._seen    = set()    # figure paths already added
        self._subtitle = ''

    # ------------------------------------------------------------------ #
    # Builders
    # ------------------------------------------------------------------ #
    def set_subtitle(self, text):
        self._subtitle = str(text)
        return self

    def add_meta(self, label, value):
        self._meta.append((str(label), _fmt(value)))
        return self

    def add_metric(self, label, value, unit=''):
        self._metrics.append((str(label), _fmt(value), str(unit)))
        return self

    def add_section(self, title):
        self._items.append({'type': 'section', 'title': str(title)})
        return self

    def add_table(self, title, headers, rows, note=''):
        '''Add a table; ``rows`` is a sequence of sequences.'''
        body = []
        for row in rows:
            body.append('<tr>%s</tr>' % ''.join(
                '<td>%s</td>' % _fmt(v) for v in row))
        table = ('<table><thead><tr>%s</tr></thead><tbody>%s</tbody></table>'
                 % (''.join('<th>%s</th>' % html.escape(str(h)) for h in headers),
                    ''.join(body)))
        self._items.append({'type': 'table', 'title': str(title),
                            'html': table, 'note': str(note)})
        return self

    def add_figure(self, path, caption=''):
        '''Add one image; ``path`` may be absolute or relative to outdir.'''
        full = path if os.path.isabs(path) else os.path.join(self.outdir, path)
        full = os.path.abspath(full)
        if full in self._seen or not os.path.isfile(full):
            return self
        self._seen.add(full)
        cap = caption or os.path.splitext(os.path.basename(full))[0]
        self._items.append({'type': 'figure', 'path': full, 'caption': str(cap)})
        return self

    def add_download(self, path, label=None):
        full = path if os.path.isabs(path) else os.path.join(self.outdir, path)
        if not os.path.isfile(full):
            return self
        self._items.append({'type': 'download', 'path': os.path.abspath(full),
                            'label': str(label or os.path.basename(full))})
        return self

    def add_config(self, text, title='Input configuration'):
        if not text:
            return self
        self._items.append({'type': 'config', 'title': str(title),
                            'text': str(text)})
        return self

    def add_figures_from_dir(self, patterns=FIGURE_PATTERNS, directory=None):
        '''Add every known figure found in ``directory`` (default: outdir).'''
        directory = os.path.abspath(directory or self.outdir)
        for caption, pattern in patterns:
            matches = sorted(glob.glob(os.path.join(directory, pattern)))
            matches = [m for m in matches if m.lower().endswith(
                ('.png', '.jpg', '.jpeg', '.gif', '.svg', '.webp'))]
            for m in matches:
                stem = os.path.splitext(os.path.basename(m))[0]
                text = caption if len(matches) == 1 else '%s (%s)' % (caption, stem)
                self.add_figure(m, text)
        return self

    def add_downloads_from_dir(self, patterns=DOWNLOAD_PATTERNS, directory=None):
        '''Add result files found in ``directory`` (default: outdir).'''
        directory = os.path.abspath(directory or self.outdir)
        for pattern in patterns:
            for m in sorted(glob.glob(os.path.join(directory, pattern))):
                if os.path.basename(m) == 'report.html':
                    continue
                self.add_download(m)
        return self

    # ------------------------------------------------------------------ #
    # Rendering
    # ------------------------------------------------------------------ #
    def _href(self, path):
        try:
            return os.path.relpath(path, self.outdir).replace(os.sep, '/')
        except Exception:
            return os.path.basename(path)

    def _render_figure(self, item):
        path, caption = item['path'], html.escape(item['caption'])
        name = os.path.basename(path)
        if self.embed:
            ext = os.path.splitext(path)[1].lower()
            mime = _IMAGE_MIME.get(ext, 'application/octet-stream')
            with open(path, 'rb') as fid:
                data = base64.b64encode(fid.read()).decode('ascii')
            img = 'data:%s;base64,%s' % (mime, data)
            return ('<figure><img src="%s" alt="%s"><figcaption>%s</figcaption>'
                    '</figure>' % (img, caption, caption))
        href = self._href(path)
        return ('<figure><a href="%s"><img src="%s" alt="%s"></a>'
                '<figcaption>%s</figcaption></figure>'
                % (href, href, caption, caption))

    def _render_downloads(self, items):
        rows = []
        for it in items:
            path = it['path']
            size = os.path.getsize(path) / 1024.0
            unit = 'kB' if size < 1024 else 'MB'
            val = size if size < 1024 else size / 1024.0
            rows.append('<li><a href="%s">%s</a><span>%.1f %s</span></li>'
                        % (self._href(path), html.escape(it['label']), val, unit))
        return '<ul class="dl">%s</ul>' % ''.join(rows)

    def _render_body(self):
        out = []
        if self._metrics:
            cards = ''.join(
                '<div class="card"><div class="v">%s<span class="u">%s</span></div>'
                '<div class="k">%s</div></div>'
                % (v, (' ' + u) if u else '', html.escape(k))
                for k, v, u in self._metrics)
            out.append('<div class="cards">%s</div>' % cards)

        pending_dl = []
        i = 0
        while i < len(self._items):
            it = self._items[i]
            kind = it['type']
            if kind == 'download':
                pending_dl.append(it)
                i += 1
                continue
            if pending_dl:
                out.append('<h2>Downloads</h2>%s' % self._render_downloads(pending_dl))
                pending_dl = []
            if kind == 'section':
                out.append('<h2>%s</h2>' % html.escape(it['title']))
            elif kind == 'table':
                out.append('<h2>%s</h2>%s%s'
                           % (html.escape(it['title']), it['html'],
                              ('<p class="note">%s</p>' % html.escape(it['note']))
                              if it['note'] else ''))
            elif kind == 'figure':
                # collect consecutive figures into one grid
                figs = [it]
                j = i + 1
                while j < len(self._items) and self._items[j]['type'] == 'figure':
                    figs.append(self._items[j])
                    j += 1
                out.append('<div class="figs">%s</div>'
                           % ''.join(self._render_figure(f) for f in figs))
                i = j
                continue
            elif kind == 'config':
                out.append('<h2>%s</h2><details><summary>Show</summary>'
                           '<pre>%s</pre></details>'
                           % (html.escape(it['title']), html.escape(it['text'])))
            i += 1
        if pending_dl:
            out.append('<h2>Downloads</h2>%s' % self._render_downloads(pending_dl))
        return '\n'.join(out)

    def _render_meta(self):
        if not self._meta:
            return ''
        return '<div class="meta">%s</div>' % ''.join(
            '<span>%s: <b>%s</b></span>' % (html.escape(k), v)
            for k, v in self._meta)

    def write(self, filename='report.html'):
        '''Render the report and return its path.'''
        from string import Template
        os.makedirs(self.outdir, exist_ok=True)
        page = Template(_load_template()).safe_substitute(
            title=html.escape(self.title),
            subtitle=html.escape(self._subtitle),
            meta=self._render_meta(),
            body=self._render_body())
        path = os.path.join(self.outdir, filename)
        with open(path, 'w', encoding='utf-8') as fid:
            fid.write(page)
        logging.info('HTML report written to %s.', path)
        return path
