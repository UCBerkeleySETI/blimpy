## Using `blimpy` with `docker`

### Pulling an image

```bash
docker pull <repo>:<tag>
```
Avaiable tags:
- py3_manual_stable
- py3_kern_stable
- py2_manual_stable
- py2_kern_stable

### Run container from image

```bash
docker run --name <name> -it <repo>:<tag> bash
```
This command opens an interactive bash session inside a container and assigns its name to `<name>`. You can basically treat the container as a remote machine.

### Misc Commands

To disconnect from the container, press: ctrl+P then ctrl+Q
To reconnect: `docker attach <name>`
Copying files to container: `docker cp local_file <name>:path_inside_container`
Copying files from container: `docker cp <name>:file_in_container local_path`
