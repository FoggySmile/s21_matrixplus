FROM ubuntu

RUN apt-get update
RUN apt-get install -y gcc
RUN apt-get install -y valgrind
RUN apt-get install -y make
RUN apt-get install -y vim
RUN apt-get install -y check
RUN apt-get install -y clang-format
RUN apt-get install -y lcov
RUN apt-get -y install g++
RUN apt-get install -y g++ libgtest-dev


COPY . /project
