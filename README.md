# Cellranger
This repository provides a workflow for running Cellranger, a bioinformatics tool for single-cell gene expression analysis. Below are instructions for building and running the workflow using various tools and platforms.

- Build and Push Docker Image
To build and push the Cellranger Docker image, execute the following command:

```shell
sh docker_build.sh
```
This command runs the docker_build.sh script, which builds the Docker image and pushes it to the specified registry.

- Validating WDL File
To validate the Cellranger WDL file, use the following command:

```shell
java -jar bin/womtool-59.jar validate cellranger.wdl
```
This command validates the syntax and validity of the Cellranger WDL file using the womtool tool.

- Generating input.json for Workflow
To generate the input.json file for the workflow, run the following command:

```shell
java -jar bin/womtool-59.jar inputs cellranger.wdl > input.json
```
This command generates the input template in JSON format for the Cellranger workflow using the womtool tool.

- Running WDL Locally on Cromwell Engine
To execute the Cellranger workflow locally using the Cromwell engine, use the following command:

```shell
java -jar bin/cromwell-58.jar run cellranger.wdl --inputs input.json
```
This command runs the Cellranger workflow locally with the specified input JSON file using the Cromwell engine.


## Related Links
- WDL Documentation: https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md
- Dockstore: https://dockstore.org/
- Docker Hub: https://docs.docker.com/docker-hub/
- GATK Workflows: https://github.com/gatk-workflows/
- Learn WDL: https://github.com/openwdl/learn-wdl\

Please feel free to modify the content further and add any additional information that is relevant to your project.

