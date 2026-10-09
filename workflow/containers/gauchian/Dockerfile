FROM python:3.11.11-slim-bookworm

# PyPI gauchian 1.0.2 is the pin. No public biocontainer tag was available.
# The Gauchian task pip-installs this same version when the image has no
# gauchian binary. Push this image and set gauchian_docker to skip that install.
RUN pip install --no-cache-dir gauchian==1.0.2
