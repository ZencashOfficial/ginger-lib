FROM ubuntu:18.04

MAINTAINER infrastructure@zensystem.io

SHELL ["/bin/bash", "-c"]

COPY ./ci/* /usr/local/bin/
COPY . /ginger-lib
ENV CONTAINER_JAVA_VER=openjdk-11-jdk
ENV CONTAINER_RUST_VER=nightly

# Get Ubuntu packages
RUN set -eux && export GOSU_VERSION=1.12 && export DEBIAN_FRONTEND=noninteractive \
    && apt-get update \
    && apt-get install -y --no-install-recommends build-essential ca-certificates curl dirmngr "$CONTAINER_JAVA_VER" \
       gcc-mingw-w64-x86-64 gnupg2 wget; \
# save list of currently installed packages for later so we can clean up
    savedAptMark="$(apt-mark showmanual)"; \
    apt-get update; \
    apt-get install -y --no-install-recommends ca-certificates wget ; \
    if ! command -v gpg; then \
      apt-get install -y --no-install-recommends gnupg2 dirmngr; \
    elif gpg --version | grep -q '^gpg (GnuPG) 1\.'; then \
# "This package provides support for HKPS keyservers." (GnuPG 1.x only)
      apt-get install -y --no-install-recommends gnupg-curl; \
    fi; \
    rm -rf /var/lib/apt/lists/*; \
    \
    dpkgArch="$(dpkg --print-architecture | awk -F- '{ print $NF }')"; \
    wget -O /usr/local/bin/gosu "https://github.com/tianon/gosu/releases/download/$GOSU_VERSION/gosu-$dpkgArch"; \
    wget -O /usr/local/bin/gosu.asc "https://github.com/tianon/gosu/releases/download/$GOSU_VERSION/gosu-$dpkgArch.asc"; \
    \
# clean up fetch dependencies
    apt-mark auto '.*' > /dev/null; \
    [ -z "$savedAptMark" ] || apt-mark manual $savedAptMark; \
    apt-get purge -y --auto-remove -o APT::AutoRemove::RecommendsImportant=false; \
    \
    chmod +x /usr/local/bin/gosu; \
# verify that the binary works
    gosu --version; \
    gosu nobody true \
    && chmod +x /usr/local/bin/entrypoint.sh \
    && chmod +x /usr/local/bin/build.sh \
    && apt-get -y clean \
    && apt-get -y autoclean \
    && rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*.deb

RUN set -eux \
# install rust
    && curl https://sh.rustup.rs -sSf | sh -s -- -y --default-toolchain nightly \
    && source /root/.cargo/env \
    && rustup default $CONTAINER_RUST_VER \
    && rustup toolchain install "$CONTAINER_RUST_VER" \
    && rustup target add --toolchain "$CONTAINER_RUST_VER" x86_64-pc-windows-gnu \
    && mkdir -p  $HOME/.cargo \
    && echo "[target.x86_64-pc-windows-gnu]" > $HOME/.cargo/config \
    && echo 'linker = "x86_64-w64-mingw32-gcc"' >> $HOME/.cargo/config \
    && echo 'source $HOME/.cargo/env' >> $HOME/.bashrc

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

