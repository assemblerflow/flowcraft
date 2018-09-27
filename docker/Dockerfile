FROM python:3.6-alpine3.7
MAINTAINER Bruno Gonçalves <bfgoncalves@medicina.ulisboa.pt>

RUN apk add --no-cache git

WORKDIR /flowcraft

# Clone FlowCraft
RUN git clone https://github.com/assemblerflow/flowcraft.git
WORKDIR ./flowcraft

# Install flowcraft
RUN python setup.py install

WORKDIR /flowcraft

# Remove unnecessary packages
RUN apk del git