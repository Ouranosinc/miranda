"""Validate URLs."""

from __future__ import annotations
import re
import typing


__all__ = ["url_validate"]


def url_validate(target: str) -> typing.Match[str] | None:
    """
    Validate whether a supplied URL is reliably written.

    Parameters
    ----------
    target : str
        The URL to validate.

    Returns
    -------
    typing.Match[str], optional
        The match object if the URL is valid.

    References
    ----------
    https://stackoverflow.com/a/7160778/7322852
    """
    url_regex = re.compile(
        r"^(?:http|ftp)s?://"  # http:// or https://
        # domain...
        r"(?:(?:[A-Z\d](?:[A-Z\d-]{0,61}[A-Z\d])?\.)+(?:[A-Z]{2,6}\.?|[A-Z\d-]{2,}\.?)|"
        r"localhost|"  # localhost...
        r"\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})"  # ...or ip
        r"(?::\d+)?"  # optional port
        r"(?:/?|[/?]\S+)$",
        re.IGNORECASE,
    )
    return re.match(url_regex, target)
