FROM continuumio/miniconda3

RUN conda update -q -y conda
RUN conda install -q -y wget
RUN apt-get install -y procps
RUN conda install -c conda-forge/label/gcc7 r-htmlwidgets
RUN conda install -c conda-forge/label/gcc7 r-networkd3
RUN conda install -c conda-forge/label/gcc7 r-tibble
RUN conda install -c conda-forge/label/cf201901 pandoc
RUN conda install -c conda-forge/label/gcc7 r-webshot
#RUN conda install -c conda-forge/label/cf201901 phantomjs
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
      ca-certificates \
      bzip2 \
      libfontconfig \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
      curl \
  && mkdir /tmp/phantomjs \
  && curl -L https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-linux-x86_64.tar.bz2 \
          | tar -xj --strip-components=1 -C /tmp/phantomjs \
  && cd /tmp/phantomjs \
  && mv bin/phantomjs /usr/local/bin \
  && cd \
  && apt-get purge --auto-remove -y \
      curl \
  && apt-get clean \
  && rm -rf /tmp/* /var/lib/apt/lists/*

# Run as non-root user
RUN useradd --system --uid 72379 -m --shell /usr/sbin/nologin phantomjs

USER phantomjs
EXPOSE 8910
CMD ["/usr/local/bin/phantomjs"]
# PhantomJS requires  an OpenSSL config even if it's an empty one,
# else it'll complain about "libssl_conf.so: cannot open shared object file"
# which seems to be a recent bug.
ENV OPENSSL_CONF=/opt/openssl.cnf