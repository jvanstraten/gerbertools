# Manylinux version to use.
ARG MANYLINUX=2014
FROM quay.io/pypa/manylinux${MANYLINUX}_x86_64

# Max number of parallel processes to use.
ARG NPROCS=100
ENV NPROCS ${NPROCS}

# Install auditwheel.
ARG PYTHON_VERSION=38
ENV PYBIN /opt/python/cp${PYTHON_VERSION}-cp${PYTHON_VERSION}*/bin
RUN ${PYBIN}/pip3 install -U pip auditwheel

WORKDIR /src

ADD . .
