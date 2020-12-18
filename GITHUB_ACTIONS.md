This repository uses GitHub Actions workflows to test, validate, and publish code. The configuration files are listed inside `.github/workflows`.

### Workflow

#### On Commit or Pull-Request
Test and validate the integrity of each commit to any branch.

1. `python_tests.yml`: Run Python tests with coverage report.
2.  `docker_build.yml`: Run build test with Docker.

#### On Master Commit
Publish the image to Docker Hub after a commit to `master` branch.

1. `push_docker.yml`: Build & publish the image on Docker Hub.

### Required Secrets
- **DOCKER_USER**: Docker Hub Username.
- **DOCKER_PASS**: Docker Hub Password.
- **CODECOV_TOKEN**: Codecov Blimpy Token.
