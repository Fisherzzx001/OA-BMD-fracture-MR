#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
317_ssh_proxy_relay.py -- minimal HTTP-CONNECT relay for OpenSSH ProxyCommand.

The sandbox permits CONNECT to ssh.github.com:443 but blocks git's HTTPS
smart-HTTP on github.com:443.  This helper lets `git push` use GitHub SSH
(ssh.github.com:443) through the sandbox HTTP proxy.

Usage (from the OA_B_PROJECT root, via GIT_SSH_COMMAND):

  ProxyCommand='python3 code/317_ssh_proxy_relay.py ssh.github.com 443'

The proxy host/port is taken from HTTPS_PROXY (falling back to
127.0.0.1:49831).  stdin/stdout are relayed verbatim; nothing is logged, so
SSH credentials never appear in output.
"""
import os
import re
import sys
import socket
import threading

DEFAULT_PROXY = "127.0.0.1:49831"


def proxy_addr():
    u = os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy") or ""
    m = re.search(r"(?:https?://)?([^:/@]+):(\d+)", u)
    if m:
        return m.group(1), int(m.group(2))
    host, port = DEFAULT_PROXY.split(":")
    return host, int(port)


def relay(sock):
    """Bidirectional copy between the tunnel socket and stdin/stdout."""
    inp = sys.stdin.buffer
    out = sys.stdout.buffer

    def c2s():
        try:
            while True:
                # read1() returns as soon as ANY bytes are available; a plain
                # read(n) would block until n bytes arrive and stall the KEX.
                b = inp.read1(65536)
                if not b:
                    break
                sock.sendall(b)
        except (OSError, ValueError):
            pass
        # Deliberately do NOT half-close the socket on stdin EOF: GitHub's SSH
        # server sends its banner first, and an immediate SHUT_WR can make the
        # proxy tear the tunnel down before the handshake starts.  OpenSSH
        # terminates this ProxyCommand process when the session ends anyway.

    def s2c():
        try:
            while True:
                b = sock.recv(65536)
                if not b:
                    break
                out.write(b)
                out.flush()
        except (OSError, ValueError):
            pass

    t = threading.Thread(target=c2s, daemon=True)
    t.start()
    s2c()
    t.join(timeout=2)


def main():
    host = sys.argv[1] if len(sys.argv) > 1 else "ssh.github.com"
    port = int(sys.argv[2]) if len(sys.argv) > 2 else 443
    ph, pp = proxy_addr()
    s = socket.create_connection((ph, pp), timeout=20)
    s.settimeout(60)
    req = (f"CONNECT {host}:{port} HTTP/1.1\r\n"
           f"Host: {host}:{port}\r\n"
           f"Proxy-Connection: Keep-Alive\r\n\r\n")
    s.sendall(req.encode("ascii"))
    buf = b""
    while b"\r\n\r\n" not in buf:
        chunk = s.recv(4096)
        if not chunk:
            raise RuntimeError("proxy closed during CONNECT")
        buf += chunk
    head = buf.split(b"\r\n\r\n", 1)[0].decode("latin-1", "replace")
    if not re.search(r"HTTP/1\.[01] 200", head):
        raise RuntimeError("proxy CONNECT failed: " + head.split("\r\n")[0])
    # Any bytes after the header are already SSH server data; push them out.
    extra = buf.split(b"\r\n\r\n", 1)[1]
    if extra:
        sys.stdout.buffer.write(extra)
        sys.stdout.buffer.flush()
    s.settimeout(None)
    relay(s)


if __name__ == "__main__":
    main()
