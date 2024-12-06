
This document provides a guide to replicate the experiments described in the paper:

> *Present/Absent $k$-mer Alignment-Free Dissimilarity Measures  for Genome Scale Analysis: Statistical Power, False Positives Control  and Combinatorial Limits*

We detail the steps necessary to prepare the environment, generate datasets, compute alignment-free (AF) measures, and analyze results. 

# Initial Requirements

As a requirement, we are assuming the availability of:
- a _Java_ compliant virtual machine (>= 1.8, [https://adoptopenjdk.net](https://adoptopenjdk.net/)).
- the _Spark_ framework (>= 3.5.1,[https://spark.apache.org/downloads.html](https://spark.apache.org/downloads.html) ) 
- a _Python_ interpreter (>= 3.8, [https://python.org/](https://python.org/))
- the _MAVEN_ framework (>=3.3, [https://maven.apache.org/download.cgi](https://maven.apache.org/download.cgi))
- the _R_ framework (>=4.1, [https://www.r-project.org/](https://www.r-project.org/))
- the _KMC_ k-mer counting tool (>=3, [https://github.com/refresh-bio/KMC](https://github.com/refresh-bio/KMC))
- the _Mash_ distance estimator tool (>=2.3,[https://github.com/marbl/Mash])


Once downloaded and unpacked the archive file, a new directory will be created with name **ALPACA-main**.  In order to build all components needed to execute the experiments, move in this directory and run the command:

> mvn package

As a result, MAVEN will start downloading all required libraries and building the code needed to run the different steps of the experiment. At the end of the process, a new **ALPACA-1.0-SNAPSHOT.jar** file will be available in the **target** directory.



# Usage

The **ALPACA** package provides a set of tools for executing the different steps required to evaluate the control of Type I error and the power of the test statistic over a set of AF functions under investigation, by means of a collection of synthetic datasets generated according to different null and alternate models. In the following, we report the instructions needed to carry out the different steps of this procedure.

## Step 1: Dataset generation

In this step, new datasets are generated for the study of histogram-based AF functions. 

It is possible to generate the same datasets considered in the paper, as well as other ones, by using the dataset generator available in this package. It is able to generate pair of sequences of increasing size according to different generative models. This is done according to the following procedure:



 1. Given an alphabet $\Phi=\{A,T,G,C\}$ with probability  $p=\{0.25,0.25,0.25,0.25\}$  and two integers  $n$ and $m$, $m$ pairs of sequences  each of length $n$ are generated using random samplings from a multinomial distribution with parameters $p$ and $n$, denoted  $Multinom(p,n)$. We refer to this model as $NM$.
 
 2. Fix an integer $\ell$, whose role will be clear immediately, each pair  of sequences $(x,y)$, generated in the previous stage, is made more similar as follows.  For each position $i$ of $x$ (with $i=1\ldots n$), we randomly sample from a Bernoulli distribution with parameter of success $\gamma$, $Bern(\gamma)$: in case of success at position $i$, the sub-sequence $x[i, i+\ell-1]$ replaces the sub-sequence $y[i, i+\ell-1]$ and the next positions is resumed to position $i+\ell$ (the sequence that has been replaced is hereafter called motif and $\ell$ denotes its length). The result is a pair of sequences $(x,\hat{y})$, where $\hat{y}$ is intuitively more similar than $y$ to $x$. Such a similarity increases with the increase of the probability of success $\gamma$. Finally, the set of the new $m$ pairs of $(x,\hat{y})$ is returned as output. We refer to this model as $PT(\ell, \gamma)$.
 
 3. Each pair  of sequences $(x,y)$, generated in stage 1, is made more similar  by implanting a  motif in both sequences in the  position selected by  $Bern(\gamma)$.  First, $m$ pairs of sequences of length $n$ are generated as in the case of $PT$. Second, we build a set $M$ of $d$ distinct motifs, each of length $\ell$ using a simple variant of $NM$ (details left to the reader). The number $d$ is chosen to respect the proportion used in \cite{2009reinertalignment} regarding  the number of motifs available for  insertion into a pair of sequences, with respect to their length.  In that case, a single motif of length $5$ for sequences of length up to $\hat{m}=25000$ was used. Accordingly, taking as reference the maximum sequence length used in Reinert et al.,  we choose $d=\frac{n}{\hat{m}\ell}$. Third, when a position $i$ in a pair of random sequences $(x,y)$ is chosen for replacement (according to the $Bern(\gamma)$ distribution), a motif $t$ is randomly sampled from $M$ (with replacement) and it is copied in both $x[i, i+\ell-1]$ and  $y[i, i+\ell-1]$. 
We refer to this model as $MR(d, \ell, \gamma)$.

Note: the aforementioned procedure is repeated several times by considering increasing values of $n$, starting from a $base$ value and up to a $destination$ value, by means of a user-defined $stepping size$. 


The dataset generator can be recalled by executing the **it.unisa.di.bio.DatasetBuilder** class available in the **ALPACA-1.0-SNAPSHOT.jar** package.  Moreover, being a Spark application, it has to be executed through the **spark-submit** command, using the following syntax: 


    spark-submit  --class it.unisa.di.bio.DatasetBuilder target/ALPACA-1.0-SNAPSHOT.jar output_directory dataset_type [local|yarn] initial_length max_lenght step_size

where:

1. **`--output <output_directory>`**:    
   enotes the file path where to save the generated sequences. It can be either a local file system path (if *local* execution mode is chosen for Spark) or an HDFS file system path (if  *yarn* execution mode is chosen for Spark)  

2. **`--generator <dataset_type>`**:    
   Indicates the type of dataset to generate.  
   - Currently, only **`detailed`** is supported.  

3. **`--mode <execution_mode>`**:    
   Specifies whether the code should be run locally or on an Hadoop cluster.  
   - Accepted values are **`local`** or **`yarn`**.  

4. **`--from-len <initial_length>`**:   
   Denotes the initial length (in characters) of the sequences to be generated.

5. **`--to-len <max_length>`**:   
   Specifies the maximum length (in characters) of the sequences to be generated.

6. **`--step <step_size>`**: 
   Denotes the length increase, in characters, of all sequences to be generated  with respect to the previous generation

7. **`--pairs <numPairs>`**:  
   Specifies the number of sequence pairs to generate for each length.

8. **`--geneSize <geneSize>`**:  
   Defines the size of the gene fragments to shuffle.

9. **`--patternSize <patternSize>`**:  
   Indicates the length of the motif or pattern used in the generation process.

A convenience script, called **runDatasetBuilder.sh**, is available in the  **scripts** directory for simplifying its execution. It accepts eight input parameters that must be provided during execution. Below is a detailed explanation of each parameter:


1. **`model`** (required):  
   Specifies the generative model to use for creating sequences. Supported values are:  
   - **`uniform`**: Uses the uniform random model (`NM`) to generate sequences.  
   - **`eColiShuffled`**: Uses a model based on shuffled *E. coli* genome sequences.  

2. **`fromLen`** (required):  
   Denotes the starting sequence length (in characters) for dataset generation.

3. **`toLen`** (required):  
   Specifies the maximum sequence length (in characters) to be generated.

4. **`step`** (required):  
   Denotes the length increase, in characters, of all sequences to be generated  with respect to the previous generation

5. **`#pairs`** (required):  
   Specifies the number of sequence pairs to generate for each length.

6. **`geneSize`** (required):  
   Defines the size of the gene fragments to shuffle. 

7. **`patternSize`** (required):  
   Indicates the length of the pattern used in the generation process.

8. **`executionMode`** (required):  
   Specifies the execution environment for Spark. Supported values are:  
   - **`local`**: Runs the generator locally.  
   - **`yarn`**: Executes the generator on a Hadoop cluster.

We advise to modify and use this script, according to the configuration of the Spark cluster being used. In this same script, it is required to provide the name of the existing directory where to save the generated sequences. By default, this script is assumed to be run from the package main directory.

---

### Example Usage

```bash
./runDatasetBuilder.sh uniform 1000 10000 2000 1000 1 32 local
```

Detailed instructions about the datasets to generate are provided by means of a configuration file called **PowerStatistics.properties**, that has to be found in the local directory where the generation script is executed. Details about the options supported by the dataset generator, together with a live example useful to recreate the same type of datasets considered in the paper, is provided in the **PowerStatistics.properties.sample** file available in the **src/main/resources** package directory. We advise to make a copy of this file, to be renamed in 
**PowerStatistics.properties**, and modify it according to the datasets to be generated.

If run with the sample **PowerStatistics.properties** file, the generator will return as output four different collection of files: 

    Uniform.*: FASTA files generated using the $NM$ model. The file names will reflect the number of pairs of sequences, the maximum length of each sequence.
    Uniform-T1.*: FASTA files generated using the $NM$ model for Type I error control.  The file names will reflect the number of pairs of sequences, the maximum length of each sequence.
    PatTransf-U-*: FASTA files generated using the $PT$ model. The file names will reflect the number of pairs of sequences, the maximum length of each sequence and the value used for the parameter \gamma.
    MotifRepl-U-*: FASTA files generated using the $MR$ model. The file names will reflect the number of pairs of sequences, the maximum length of each sequence and the value used for the parameter \gamma.

# Step 2: AF present/absent measures evaluation
Let $S$ be a set of $m$ pairs of sequences, each of length $n$,  generated in the previous step via $NM$. Let $S^{PT}\_{\gamma}$ be a set of $m$ pairs of sequences, generated starting from $S$ during Step 1,  via $PT$ with parameter $\gamma$. Let $S^{MR}\_{\gamma,m}$ be a set of $m$ pairs of sequences, generated starting from $S$ during Step 1,  via $MR$ with parameter $\gamma$ and $m$.

In this step, we evaluate several types of AF  measures over all distinct pairs of sequences belonging to $S$, $S^{PT}\_{\gamma}$  and $S^{MR}\_{\gamma,m}$, using different assignments for $k$.  
It is possible to evaluate all AF measures for all sequence pairs stored in the **target directory** used in the previous step, by executing the **PySparkPresentAbsent4.py** script.

A convenience script, called **runPresentAbsent.sh**, is available in the  **scripts** directory for simplifying its execution.  Below is a description of its parameters:

1. **`seqLen`**:  
   Specifies the sequence length (in characters) to be analyzed. This should match one of the sequence lengths generated in the previous step.

2. **`mainDataDir`**:  
   The directory where the generated datasets are stored.  
   - This should point to the root directory containing the sequence data to be analyzed.

3. **`executionMode`**:  
   Specifies the Spark execution environment.  
   - Accepted values are:
     - **`local`**: Executes the script on the local machine.
     - **`yarn`**: Executes the script on a YARN-based Hadoop cluster.  

   
Notice that if `yarn` mode is selected, additional configuration parameters will be used for instrumenting the Hadoop execution. Further information about these parameters can be found in the **runPresentAbsent.sh** file.
Additionally, the script supports parameters for managing external tools such as `Mash` and `KMC`, allowing users to specify the paths where these tools are installed. Users can also define specific directories for input datasets and output files, providing control over the organization and accessibility of generated results. 

# Step 3: AF present/absent measures summarization

In this step, we first clean and reorganize the output generated in the previous step by producing a single  CSV file containing all AF present/absent measures evaluated during the process.

This is achieved by running either the **`makeCSVReport.sh`** or **`makeCSVReportLocal.sh`** script, depending on the execution mode used in the previous steps:

- Use **`makeCSVReportLocal.sh`** if the previous steps were executed in local mode.
- Use **`makeCSVReport.sh`** if the previous steps were executed in yarn mode (i.e., using the HDFS file system)


For the `makeCSVReportLocal.sh` script, use the following syntax:

```bash
./makeCSVReportLocal.sh inputDirectory outputFile
``` 
Below is a description of its parameters:

1. **`inputDirectory`**:  
   The path to the folder containing the CSV files generated in the previous step
2. **`outputFile`**:  
   The path of the CSV file that will be generated as output  


For the `makeCSVReport.sh` script, use the following syntax. Note that this script assumes the CSV files
generated during the previous step have already been transferred to the local file system before execution.

```bash
./makeCSVReport.sh inputDirectory outputFile
```
Below is a description of its parameters:

1. **`inputDirectory`**:  
   The path to the folder containing the CSV files generated in the previous step
2. **`outputFile`**:  
   The path of the CSV file that will be generated as output  




evaluate the control of Type I error and the power of the test statistic over a set of input AF functions by analyzing the AF values determined in the previous step respectively on the PT and MR alternate models, and on the null distribution, with different values of $k$, $\alpha$, and $\gamma$. Then, we summarize and plot the results of this evaluation on several charts that are automatically saved on disks as graphical images. The whole analysis process takes place in three substeps.

## Powerstatistics evaluation
In this substep, we evaluate the power statistics of choice on the AF measures determined in the previous step. This is done by running the **it.unisa.di.bio.powerstatistics.PowerEvaluator** class available in the **powerstatistics-1.0-SNAPSHOT.jar** package, using the following syntax:

    java -cp ALPACA-1.0-SNAPSHOT.jar  it.unisa.di.bio.powerstatistics.PowerEvaluator input_directory [MotifReplace|PatternTransfer|Both] syntheticAllLen


Here, the first argument reports the directory containing the AF measures evaluated during previous step, encoded as CSV files. The second argument defines which alternate model to consider. The third argument must be set to *syntheticAllLen*.
As a result, the program will generate, for each AM,  each value of $\alpha$ and of $k$, a separate directory containing a set of JSON files containing the corresponding power statistics.

## Powerstatistics summarization 
In this substep, we summarize the results available in the JSON files produced in the previous substep, as well as the CSV files containing the AF measures evaluated in the previous step, in a pair of **R** dataframes ready for the analysis.  This is done by running the following command lines, to be executed from the package home directory:

    R -f R-Scripts/Power+T1-Json2RDS.R
    R -f R-Scripts/RawDistances-CSV2RDS.R

Once run, these two scripts will produce two files, **RawDistances-All.RDS** and **Power+T1-Results.RDS**, used as reference during the next substep related to charts' generation. 

## Powerstatistics charting
In this substep, data available in files  **RawDistances-All.RDS** and **Power+T1-Results.RDS** is processed and presented by means of several types of charts. In the following we report the list of plotting scripts available with a short description about their expected output:

[//]: # ( - **Plot-T1OnlyPanel.R**:  Shows several different collections of colored boxplots &#40;one for ech distinct value of $k$&#41; reporting the results of the T1 Error Control experiments for each of the considered AF functions &#40;on the abscissa&#41;. The vertical length of each boxplot is proportional to the percentage of false positives. )

[//]: # ( Used to generate **Figure 1** of the main paper, and **Figure 1** of the supplementary material.)

[//]: # ()
[//]: # ( - **PlotHammingDistances.R**: Shows a collection of colored boxplots reporting the Hamming Distance of values generated by the Alternative models, with respect to the Null model. On the abscissa, it is shown the generative models considered, with the letter $G$ denoting $\gamma$.)

[//]: # (Used to generate **Figure 2** of the main paper. )

[//]: # ()
[//]: # ( -  **Heatmap-Cluster.R**: Creates a heatmap of the delta values obtained as the difference between the average of the distribution of AF values computed with   $NM$ and one of $MR$ and $PT$, with for different combinations of $n$, $k$ and $\gamma$,  reported as colors on the left panel annotation. The dendrogram on the top is a hierarchical clustering  on the delta values with Euclidean distance and complete linkage.)

[//]: # (Used to generate **Figure 3** of the main paper. )

[//]: # ()
[//]: # (- **PlotPanelAllMeasures.R**: Produces two panels reporting the power trend, respectively for the $MR$ and $PT$ alternative models. In each panel, for each AF and $\gamma$ is reported the power level obtained across different values of $n$. It is also colored according to the value of $k$.)

[//]: # (Used to generate **Figure 4** of the main paper. )

[//]: # ()
[//]: # ()
[//]: # (- **PlotRawDistances-AllK.R**: # Produces a panel made of a set of stacked boxplots &#40;one for each distinct value of k&#41;, for each AF.)

[//]: # ( Each boxplot reports the distribution of the proportion of true positives, cumulatively by length, by considering values generated by the Alternative models, with increasing gamma, with respect to $NM$.)

[//]: # (Used to generate **Figure 5** of the main paper. )

[//]: # ()
[//]: # (- **PlotAllRawDistances.R**: As **PlotOneMeasureRawDistances-AllK.R**, but by considering multiple AF measures at same time.)

[//]: # (Used to generate **Figure 2** and **Figure 3** of the supplementary material. )

[//]: # ()
[//]: # ( - **PanelBoxplot-Power-AllMeasures.R**: Creates a collection of images, one for each considered Alternative model. In each image it is reported. for each considered AF function and each value of $\gamma$, the power levels obtained across different values of $n$, as a function of $k$.)

[//]: # (Used to generate **Figure 4** of the supplementary material. )


Each of these scripts can be executed using the following syntax:

> R -f <script name>

Notice that each of these scripts comes with a set of basic options (e.g., changing the input directory, defining the values of $\gamma$ to be considered) that are available in the **Options** section of their initial lines and that can be customized, according to the own's needs, using a text editor. 








 

