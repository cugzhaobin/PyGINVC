"""
Streamlit Logging Handler - Captures Python logging output into a list
that can be displayed in the Streamlit UI.
"""
import logging


class StreamlitLogHandler(logging.Handler):
    """
    A logging handler that appends formatted log records to a list.
    Designed to work with ``st.session_state["log_messages"]``.
    """

    def __init__(self, log_list: list):
        super().__init__()
        self.log_list = log_list
        self.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s",
                              datefmt="%H:%M:%S")
        )

    def emit(self, record):
        try:
            msg = self.format(record)
            self.log_list.append(msg)
            # Keep a reasonable limit to avoid memory issues
            if len(self.log_list) > 5000:
                self.log_list[:] = self.log_list[-5000:]
        except Exception:
            self.handleError(record)


def setup_streamlit_logging(log_list: list):
    """
    Install the StreamlitLogHandler on the root logger so that all
    pyginvc module logs are captured.

    Parameters
    ----------
    log_list : list
        A mutable list (typically ``st.session_state["log_messages"]``)
        that will receive formatted log messages.
    """
    root_logger = logging.getLogger()
    # Avoid adding duplicate handlers on Streamlit reruns
    for handler in root_logger.handlers:
        if isinstance(handler, StreamlitLogHandler):
            return
    handler = StreamlitLogHandler(log_list)
    handler.setLevel(logging.INFO)
    root_logger.addHandler(handler)
    root_logger.setLevel(logging.INFO)
