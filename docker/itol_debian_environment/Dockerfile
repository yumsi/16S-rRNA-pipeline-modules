FROM python:slim

LABEL maintainer="joshua.koopmans@radboudumc.nl"
LABEL version="1.0"
LABEL description="This image contains the zip, curl and python3 packages"

ARG USER=joshuak
ARG GROUP=unicorn

# Add a non-root user
RUN useradd -ms /bin/bash ${USER}


# Install the curl, zip and python3 packages
RUN apt-get update && apt-get install curl zip --yes --force-yes

RUN pip3 install pandas

# Set added user as default user
USER ${USER}

