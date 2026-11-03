# Per-sample assessment image.
#
# Build from the repository root, then push and set PdSample.pdvar_docker:
#   docker build -f workflow/containers/pdvar/Dockerfile -t REGISTRY/pdvar:COMMIT .
#   docker push REGISTRY/pdvar:COMMIT
#
# Without this image, leave pdvar_docker at python:3.11.11-slim-bookworm and
# set pdvar_git_commit to a commit SHA of gqworks/pdsample-terra. The Assess
# task force-installs git+https://github.com/gqworks/pdsample-terra.git@SHA.
# Do not pin Assess to "main": Cromwell call-cache keys on the input string, so
# a moving branch tip will reuse a stale Assess result. Script/config assets use
# assets_git_commit (volatile FetchRepoAssets) separately.
FROM python:3.11.11-slim-bookworm

WORKDIR /opt/src
COPY pyproject.toml README.md ./
COPY src ./src
COPY config ./config
COPY resources/panel ./resources/panel
RUN pip install --no-cache-dir .
