# concTree-NF

Phylogenetic Trees and Bootstrap Support Values Utilising Alignment Uncertainty from Multiple MSA Generated Trees

## Motivation

concTree incorporates the information from  uncertainty in multiple sequence alignments into the construction of multiple phylogenetic trees and the bootstrap support values.


## Schematic Outline



## Quick start 

Make sure you have all the required dependencies listed in the last section.

Install the Nextflow runtime by running the following command:

    $ curl -fsSL get.nextflow.io | bash


When done, you can launch the pipeline execution by entering the command shown below:

    $ nextflow run skptic/concTree-nf
    

By default the pipeline is executed against the provided example dataset. 
Check the *Pipeline parameters*  section below to see how enter your data on the program 
command line.     

## Pipeline parameters

#### `--input` 
   
* Specifies the location of the *fasta* file(s) containing sequences.
* Multiple files can be specified using the usual wildcards (*, ?), in this case make sure to surround the parameter string
  value by single quote characters (see the example below)
* By default it is set to the concTree-NF's location: `./tutorial/*.fa`

Example: 

    $ nextflow run skptic/concTree-nf --input '/home/dataset/*.fasta'

This will handle each fasta file as a seperate alignment/tree/bootstrap support



#### `--alignments` 
   
* Specifies the number of alignments to sample from the guidance2 alternative alignments.
* Min = 1, Max = 400
* By default it is set to `100`

Example: 

    $ nextflow run skptic/concTree-nf --alignments '25'


  
#### `--output`
   
* Specifies the folder where the results will be stored for the user.  
* It does not matter if the folder does not exist.
* By default is set to concTree-NF's folder: `./tutorial/results` 

Example: 

    $ nextflow run skptic/concTree-nf --output /home/user/my_results 
  


## Cluster support

concTree-NF execution relies on [Nextflow](http://www.nextflow.io) framework which provides an 
abstraction between the pipeline functional logic and the underlying processing system.

Thus it is possible to execute it on your computer or any cluster resource
manager without modifying it.

Currently the following platforms are supported:

  + Oracle/Univa/Open Grid Engine (SGE)
  + Platform LSF
  + SLURM
  + PBS/Torque


By default the pipeline is parallelized by spanning multiple threads in the machine where the script is launched.

To submit the execution to a SGE cluster create a file named `nextflow.config`, in the directory
where the pipeline is going to be launched, with the following content:

    process {
      executor='sge'
      queue='<your queue name>'
    }

In doing that, tasks will be executed through the `qsub` SGE command, and so your pipeline will behave like any
other SGE job script, with the benefit that *Nextflow* will automatically and transparently manage the tasks
synchronisation, file(s) staging/un-staging, etc.

Alternatively the same declaration can be defined in the file `$HOME/.nextflow/config`.

To lean more about the avaible settings and the configuration file read the 
[Nextflow documentation](http://www.nextflow.io/docs/latest/config.html).
  
  
Dependencies 
------------

 * Java 7+ 
 * Guidance2
 * ClustalW2
 * MAFFT
 * PRANK
 * T-Coffee
 * RAXML

