
This document provides a guide to replicate the experiments described in the paper:

> *Present/Absent $k$-mer Alignment-Free Dissimilarity Measures  for Genome Scale Analysis: Statistical Power, False Positives Control  and Combinatorial Limits*

We detail the steps necessary to prepare the environment, generate datasets, compute alignment-free (AF) measures, and analyze results. 

# Initial Requirements

As a requirement, we are assuming the availability of:

- a _Java_ compliant virtual machine (>= 1.8, [https://adoptopenjdk.net](https://adoptopenjdk.net/)).
- the _Spark_ framework (>= 3.5.1, [https://spark.apache.org/downloads.html](https://spark.apache.org/downloads.html)).
- a _Python_ interpreter (>= 3.8, [https://python.org/](https://python.org/)).
- the _MAVEN_ framework (>= 3.3, [https://maven.apache.org/download.cgi](https://maven.apache.org/download.cgi)).
- the _R_ framework (>= 4.1, [https://www.r-project.org/](https://www.r-project.org/)), including the following packages:
  - `DescTools`
  - `dplyr`
  - `ggplot2`
  - `hrbrthemes`
  - `filelock`
  - `optparse`
- the _KMC_ k-mer counting tool (>= 3, [https://github.com/refresh-bio/KMC](https://github.com/refresh-bio/KMC)).
- the _Mash_ distance estimator tool (>= 2.3, [https://github.com/marbl/Mash](https://github.com/marbl/Mash)).



Once downloaded and unpacked the archive file, a new directory will be created with name **ALPACA-main**.  In order to build all components needed to execute the experiments, move in this directory and run the command:

> mvn package

As a result, MAVEN will start downloading all required libraries and building the code needed to run the different steps of the experiment. At the end of the process, a new **ALPACA-1.0-SNAPSHOT.jar** file will be available in the **target** directory.



# Usage

## **1. Experiments on Randomly Generated Genomic Sequences**

The **ALPACA** package provides a set of tools for executing the different steps required to evaluate the control of Type I error and the power of the test statistic over a set of AF functions under investigation, by means of a collection of synthetic datasets generated according to different null and alternate models. In the following, we report the instructions needed to carry out the different steps of this procedure.

### Step 1: Dataset generation

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


The script generates its output in a subdirectory following this format:

**`data/<model>-<geneSize>-<#pairs>/len=<sequence_length>/`**


We advise to modify and use this script, according to the configuration of the Spark cluster being used. In this same script, it is required to provide the name of the existing directory where to save the generated sequences. By default, this script is assumed to be run from the package main directory.

---

#### Example Usage

```bash
./scripts/runDatasetBuilder.sh uniform 1000 10000 2000 1000 1 32 local
```

Running the aforementioned command line, the output of the script will be saved in the following path:

**`data/uniform-32-1000/len=1000/`**

Detailed instructions about the datasets to generate are provided by means of a configuration file called **PowerStatistics.properties**, that has to be found in the local directory where the generation script is executed. Details about the options supported by the dataset generator, together with a live example useful to recreate the same type of datasets considered in the paper, is provided in the **PowerStatistics.properties.sample** file available in the **src/main/resources** package directory. We advise to make a copy of this file, to be renamed in 
**PowerStatistics.properties**, and modify it according to the datasets to be generated.

If run with the sample **PowerStatistics.properties** file, the generator will return as output four different collection of files: 

    Uniform.*: FASTA files generated using the $NM$ model. The file names will reflect the number of pairs of sequences, the maximum length of each sequence.
    Uniform-T1.*: FASTA files generated using the $NM$ model for Type I error control.  The file names will reflect the number of pairs of sequences, the maximum length of each sequence.
    PatTransf-U-*: FASTA files generated using the $PT$ model. The file names will reflect the number of pairs of sequences, the maximum length of each sequence and the value used for the parameter \gamma.
    MotifRepl-U-*: FASTA files generated using the $MR$ model. The file names will reflect the number of pairs of sequences, the maximum length of each sequence and the value used for the parameter \gamma.

## Step 2: AF present/absent measures evaluation
Let $S$ be a set of $m$ pairs of sequences, each of length $n$,  generated in the previous step via $NM$. Let $S^{PT}\_{\gamma}$ be a set of $m$ pairs of sequences, generated starting from $S$ during Step 1,  via $PT$ with parameter $\gamma$. Let $S^{MR}\_{\gamma,m}$ be a set of $m$ pairs of sequences, generated starting from $S$ during Step 1,  via $MR$ with parameter $\gamma$ and $m$.

In this step, we evaluate several types of AF  measures over all distinct pairs of sequences belonging to $S$, $S^{PT}\_{\gamma}$  and $S^{MR}\_{\gamma,m}$, using different assignments for $k$.  
It is possible to evaluate all AF measures for all sequence pairs stored in the **target directory** used in the previous step, by executing the **PySparkPresentAbsent4.py** script. A convenience script, called **runPresentAbsent.sh**, is available in the  **scripts** directory for simplifying its execution.  Below is a description of its parameters:

1. **`seqLen`**:  
   Specifies the sequence length (in characters) to be analyzed. This should match one of the sequence lengths generated in the previous step.

2. **`mainDataDir`**:  
   The parent directory where the generated datasets are stored. This should point to the root directory containing the sequence data to be analyzed.

3. **`executionMode`**:  
   Specifies the Spark execution environment.  
   - Accepted values are:
     - **`local`**: Executes the script on the local machine.
     - **`yarn`**: Executes the script on a YARN-based Hadoop cluster.  

   
Notice that if `yarn` mode is selected, additional configuration parameters will be used for instrumenting the Hadoop execution. Further information about these parameters can be found in the **runPresentAbsent.sh** file.

Additionally, the **PySparkPresentAbsent4.py** script supports several other parameters required for providing the location of 
external tools such as `Mash` and `KMC`. Users can also define specific directories for input datasets and output files, as well as providing additional settings about the experiments to run, 
like the range of $k$ values to use. 


#### Example Usage

Continuing from the previous step, to evaluate AF measures for sequences of length 100 stored in the `uniform-32-4` directory, execute the following command:

```bash
./scripts/runPresentAbsent.sh 1000 data/uniform-32-1000 local
```

Running the aforementioned command, sequences in the len=100 subdirectory will be processed 
and results for AF measures evaluation will be generated in the specified environment.

After executing the script, the output is stored in the directory:

```plaintext
data/uniform-32-1000/len=1000.1000-20241212-1200.csv
```
where `1000.1000` represents the sequence length and pair count and `20241212-1200` represents the timestamp of the file creation in `YYYYMMDD-HHMM` format.
  
# Step 3: AF present/absent measures summarization and analysis
## Substep 3.1: AF present/absent measures summarization 

In this substep, we first clean and reorganize the output generated in the previous step by producing a single  CSV file containing all AF present/absent measures evaluated during the process.

This is achieved by running the **`makeCSVReport.sh`** script, using the following syntax:

```bash
./scripts/makeCSVReport.sh inputDirectory outputFile local|yarn
``` 
Below is a description of its parameters:

1. **`inputDirectory`**:  
   The path to the folder containing the CSV files generated in the previous step
2. **`outputFile`**:  
   The path of the CSV file that will be generated as output  
3. **`executionMode`**:  
   Specifies the Spark execution environment.  
   - Accepted values are:
     - **`local`**: Executes the script on the local machine.
     - **`yarn`**: Executes the script on a YARN-based Hadoop cluster.


Notice that, when run in yarn mode, this script assumes the CSV files generated during the previous step have already been transferred to the local file system before execution.

#### Example Usage

```bash
./scripts/makeCSVReport.sh data/uniform-32-1000/len=1000.1000-20241212-1200.csv data/PresentAbsentECData-uniform-32-1000.csv local
```

Assuming the output of step 2 has been saved in the **`data/uniform-32-1000/len=1000.1000-20241212-1200.csv`**, running the aforementioned command line, will cause the aggregation of all csv files existing in that directory into a single csv file called 

**`data/PresentAbsentECData-uniform-32-1000.csv`**





### Substep 3.2: Powerstatistics evaluation
In this substep, we evaluate the power statistics and the T1 error, starting from the raw data aggregated during the previous substep.
This is done by running the following command line:

```bash
Rscript R-Scripts/PresentAbsentPower+T1-CSV2RDS-parallel.R --csv path_to_input_csv --df path_to_RDS --trsh path_to_threshold_file
```

#### Parameters Description

1. **`--csv`**:  
   The input csv file containing the raw AF present/absent measures for evaluation. This file is generated during the previous step.

2. **`--df`**:  
   The output file in RDS format where the data frame containing power statistics and T1 errors will be saved.

3. **`--trsh`**:  
   The path to the csv file containing threshold values for the evaluation.

#### Example usage


```bash
Rscript R-Scripts/PresentAbsentPower+T1-CSV2RDS-parallel.R --csv data/PresentAbsentECData-uniform-32-1000.csv --df data/PresentAbsentECData-uniform-32-1000.RDS --trsh data/Threshold.csv
```

Assuming the output of step 3 has been saved in the **`data/`** directory, running the aforementioned command line, 
will imply the generation of the two files:  **`data/PresentAbsentECData-uniform-32-1000.RDS`** and  **`data/Threshold.csv`**.

### Substep 3.2: Powerstatistics Charting

In this substep, power statistics returned by previous subset is processed and presented by means of several types of charts. In the following we report the list of plotting scripts available with a short description about their expected output.
The available scripts are:

- **PresentAbsentPlot-T1+PanelPower+AN.R**:
- **PresentAbsentPlotAN.R**:

Each of these scripts can be executed using the following syntax:

```bash
Rscript R-Scripts/<script file name> --csv <input CSV file> --df <output RDS file> --dirname <output directory>
```

#### Parameters
- **`--csv <input CSV file>`**: Specifies the path to the CSV file containing the input data.  


- **`--df <output RDS file>`**: Defines the path for the RDS file where the processed data will be saved.  


- **`--outdir <output directory>`**: Specifies the directory where output plots will be saved.  



#### Example usage

To generate power statistics charts for a dataset using `PresentAbsentPlotAN.R`, execute the following command line:
```bash
Rscript R-Scripts/PresentAbsentPlotAN.R --csv data/PresentAbsentECData-uniform-32-4.csv --df data/PresentAbsentECData-uniform-32-4.RDS --outdir data
```
The script will read the input data from `data/PresentAbsentECData-uniform-32-1000.csv` and from `data/PresentAbsentECData-uniform-32-4.RDS`. 
Once run, the script generates the following charts in the `plots` directory:

- Panel-Distances-D2-AllModels.pdf
- Panel-Distances-Euclidean-AllModels.pdf
- Panel-Distances-Hamman-AllModels.pdf
- Panel-Distances-Jaccard-AllModels.pdf
- PanelADNM.pdf
- PanelAllDistancesNM.pdf
- PanelAN.pdf
- PanelANNM.pdf
- PanelD2DistancesNM.pdf
- PanelEuclideanDistancesNM.pdf

To generate power statistics charts for a dataset using `PresentAbsentPlot-T1+PanelPower+AN.R`, execute the following command line:

```bash
Rscript R-Scripts/PresentAbsentPlot-T1+PanelPower+AN.R --csv data/PresentAbsentECData-uniform-32-4.csv --df data/PresentAbsentECData-uniform-32-4.RDS --outdir data
```
The script will read the input data from `data/PresentAbsentECData-uniform-32-1000.csv` and from `data/PresentAbsentECData-uniform-32-4.RDS`. 
Once run, the script generates the following charts in the `plots` directory:

- Panel-NullModel-ANSDAAnalysis.pdf
- Panel-PatTransf-U-ANAvAnalysis.pdf
- Panel-MotifRepl-U-ANAvAnalysis.pdf
- Panel-NullModel-ANAnalysis.pdf
- PanelPowerAnalysis-PatTransf-U-A={alpha}.pdf
- PanelPowerAnalysis-MotifRepl-U-A={alpha}.pdf
- PanelPowerAnalysis-PatTransf-U-A={alpha}-G={gamma}.pdf
- PanelPowerAnalysis-MotifRepl-U-A={alpha}-G={gamma}.pdf
- T1Box-A={alpha}.pdf




## **2. Experiments on Real Genomic Sequences**

The **ALPACA** package also provides tools to evaluate alignment-free (AF) measures on real genomic sequences. 


### **Step 1: Dataset Preparation and Comparison**



In this step, real genomic sequences are processed and compared using alignment-free measures for different parameter values using the **`runSynthetic.sh`** script.
This script automates the comparison of a real genomic sequence against modified versions generated with different substitution probabilities (*theta* values). It first uses the **`makeDistance.py`** script to create additional sequences by applying random substitutions. The modified sequences are then compared to the original sequence using the **`PyPASingleSequenceOutMemory.py`** script, which evaluates alignment-free (AF) measures for specified *k*-mer values.  Results from the comparisons are logged and aggregated into a single CSV report, which includes results for all specified *theta* values and *k*-mer sizes. 


### **Usage**

The script syntax is as follows:


```bash
./runSynthetic.sh <sequence> <remoteDataDir> <local|yarn> [outputFile [theta [k]]]
```

### **Parameters**

1. **`<sequence>`** *(required)*:  
   Path to the input genomic sequence file (`.fna` file)

2. **`<remoteDataDir>`** *(required)*:  
   Directory where modified temporary outputs will be stored.

3. **`<local|yarn>`** *(required)*:  
   Specifies the execution environment:  
   - **`local`**: Runs the script on the local machine.  
   - **`yarn`**: Executes the script on a YARN-based Spark cluster.

4. **`[outputFile]`** *(optional)*:  
   Path to the final aggregated output file. If absent, a csv file with a timestamp will be generated.

5. **`[theta]`** *(optional)*:  
   A list of substitution probabilities (theta values) to introduce modifications in the sequences.   
   - Example: `"0.005 0.01 0.02 0.05"`.

6. **`[k]`** *(optional)*:  
   Specifies the *k*-mer sizes for AF analysis.
   - Example: `"4 6 8"`.

---

### **Example Usage**

To compare a *Homo sapiens* sequence locally with theta values `0.005` and `0.01` and *k*-mer sizes `4` and `6`:

```bash
./runSynthetic.sh genome1.fna Datasets local Datasets2/genome1.csv
```

The execution of the aforementioned command will imply the generation of a report, including alignment-free measures, substitution rates, and *k*-mer analysis results for each comparison.
The experiment will be run on the genome1.fna, using as k and theta values, the default ones defined in **`runSynthetic.sh`**.


### **Step 2: Powerstatistics Charting**

The **`PlotPASemiSyntheticsStatistics.R`** script generates plots for different alignment-free measures as a function of *theta* values and *k*-mer sizes.

#### **Syntax**

The script syntax is as follows:

```bash
Rscript R-Scripts/PlotPASemiSyntheticsStatistics.R --dirname <outputDirectory> --dataframe1 <distanceFile> --dataframe2 <cvFile>
```

**Parameters**:  
- **`--dirname <inputDirectory>`**: Path to the directory containing input data.  
- **`--dataframe1 <distanceFile>`**: Path to the RDS file containing distance metrics for the sequences.  
- **`--dataframe2 <cvFile>`**: Path to the RDS file containing the coefficient of variation data.
- **`--outdir <cvFile>`**: Path to the directory where the plots will be saved.



---

### **Example Usage**

In the following example, we plot the results stored in the `Datasets2` directory, about genomic sequence  `genome1`:

```bash
Rscript R-Scripts/PlotPASemiSyntheticsStatistics.R --dirname Datasets2 --outdir plots
```

Once run, the script generates the following charts in the plots directory:

   - **PanelCV-all.pdf**, plus some additional zoomed-in views
   - **PanelCVD2.pdf**, plus some additional zoomed-in views  
   - **PanelCVEuclid.pdf**, plus some additional zoomed-in views
   - **PanelD2-genome1.pdf**
   - **PanelDensities.pdf**: 
   - **PanelDensities-genome1.pdf**
   - **PanelDistances-genome1 .pdf**
   - **PanelDistances.pdf**  
   - **PanelEuclid-genome1.pdf**, plus some additional zoomed-in views  
   - **PanelRari-all.pdf**
     
---

```


 

