# Pipeline-Tools
**Pipeline-Tools** is a collection of command line tool for manipulating and summarizing VCF files.

## Installation

### Stand-alone installation
  
**Pipeline-Tools** is currently designed only for *Linux* systems. 
You will need to install and configure the following tool:  

1. [Python](https://www.python.org/) v2.7.*

    You can check your Python version by running the following command in your terminal:

    ```sh
    $ python -V
    Python 2.7.10
    ```

    To install the correct version of Python, visit the official Python [website](https://www.python.org/downloads/).

2. Python packages: *numpy*, *scipy*, *pandas*, *matplotlib*, and *[pyVCF]*

    You will need [pip](https://packaging.python.org/guides/installing-using-linux-tools/) to install the above packages.
    After installing *pip*, run the following commands in your terminal: 

    ``` sh
    # Upgrade pip
    sudo pip install -U pip
    
    pip install numpy
    pip install scipy
    pip install pandas
    pip install matplotlib
    pip install pyVCF
    ```

3. Clone the **Pipeline-Tools** repo

    ``` sh
    # clone the repo
    git clone https://github.com/alexwaldrop/Pipeline-Tools.git
    ```
    
### Docker Installation

**Pipeline-Tools** is maintained as a docker image and can be ported anywhere Docker is available.

The only pre-requisite here is the [Docker] client. Please execute the following command line to see if your system 
already have Docker-client installed or not.

```bash
$ sudo docker --version
``` 

If Docker is not installed on your system, you can get it from [Docker-client]. 

After the Docker set up, please pull the **Pipeline-Tools** Docker image from the [Docker Hub]. To do so, please run 
the following command line:

```bash
$ sudo docker pull alexwaldrop/pipeline-tools:latest
```

You can run any of the **Pipeline-Tools** modules as Docker containers as follows:

```bash
$ sudo docker run --rm --user root alexwaldrop/pipeline-tools:latest "RecodeVCF.py --help"
$ sudo docker run --rm --user root alexwaldrop/pipeline-tools:latest "CatRecodedVCF.py --help"
$ sudo docker run --rm --user root alexwaldrop/pipeline-tools:latest "SummarizeVCF.py --help"
$ sudo docker run --rm --user root alexwaldrop/pipeline-tools:latest "CatVCFSummary.py --help"

``` 

## Submodules
At present, **Pipeline-Tools** contains two tools for VCF processing:


### RecodeVCF
Transforms VCF genotype calls, annotations into searchable/sortable RecodedVCF format. 
 [Complete Documentation](./doc/RecodeVCF.md) found here.

### SummarizeVCF
Summarizes variant features 
based on annotations contained in a VCF.  [Complete Documentation](./doc/SummarizeVCF.md) found here.


## Project Status

**Pipeline-Tools** is actively under development. To request features, please contact the author listed below.

## Authors

* [Alex Waldrop](https://github.com/alexwaldrop)


[Docker]: https://www.docker.com/
[Docker-client]: https://docs.docker.com/install/
[Docker Hub]: https://hub.docker.com/
[pyVCF]:https://github.com/jamescasbon/PyVCF