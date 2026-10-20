# Docker image for build manylinux wheels for libqif

FROM quay.io/pypa/manylinux_2_28

# prefer gcc 10
ENV PATH=/usr/share/Modules/bin:/opt/rh/gcc-toolset-10/root/usr/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

ENV CMAKE_PREFIX_PATH=/opt/ortools

COPY setup.sh /tmp/
RUN bash /tmp/setup.sh