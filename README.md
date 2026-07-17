
# GRB Radio Afterglow Analysis

This project analyzes radio afterglow data from gamma ray bursts (GRBs) to 
extract spectral indices and map how GRB energetics evolve over time.

## Dependencies
Python 3.14.3 and the following libraries:
- pandas
- matplotlib
- numpy
- scipy
- astropy

Install with:
'''
pip install pandas matplotlib numpy scipy astropy
'''

## Background
This project takes data from gamma ray bursts observed in the radio band and analyzes their comparable brightness and the temporal evolution of the spectral index. Gamma ray bursts are high energy astrophysical events observed inter- and extra-galactically characterized by high energy peaked flux in gamma rays. They are divided into short(<2s) and long(>2s) bursts, of which the duration is measured by the time in which 90 percent of the peak flux is received. These events are considered to be a results of the core collapse of massive stars(for long GRBs) or from the merger of two compact objects such as neutron stars and black holes(short GRBs). The burst takes the form of a relativistic jet which is theorized due to energy constraints. Environments around GRBs are optically thick, which indicates that isotropic radiation would not be sufficient enough to produce the energetics that are observed with GRBs. Therefore, the presence of a relativistic jet satisfies energy constraints by letting the observed radiation be magnified by the squared Lorentz factor. This relativistic jet then produces an "afterglow" in the radio to X-ray bands caused by particles in the circumburst medium being swept up by the jet and subsequently radiating synchrotron. 

Analyzing the afterglow is valuable due to the ability to extract spectral indices, which measure the dependency of flux on frequency with a power law relationship of F ~ nu^α, where α is the spectral index. This relationship is assumed for synchrotron radiation because electrons follow the energy distribution 𝑑𝑁(γ) = 𝑘γ^(−𝑝)𝑑γ. Electron index(𝑝) and spectral index are related by α = −(𝑝−1)/2. In logspace, the spectral index is then measured as the slope of the flux vs frequency relationship. This slope changes at different frequencies, coined "spectral breaks". The peak frequency, minimum frequency, cooling frequency, and synchrotron self absorption frequency all contribute to calculating valuable information about the energetics of the GRB. By extracting the spectral index for distinct epochs(1.3-4.9, 4.9-8.5, 8.5-15 GHz) over time, this project can map how GRB energetics evolve over time and lead to insights about their intrinsic brightness.

## Data
The datasets from this project are of a set older GRBs dating between 1997-2003, which in their time were the brightest GRBs observed to date. I also include data from GRB130427A, a record setting burst, and GRB221009A, the brightest GRB of all time. The raw data in original .csv or .dat format is stored in data/raw. Due to the age of most of these GRBs, many datasets have invalid flux measurements and high uncertainty, which lead to insufficient spectral data. While some datasets include observations across different bands, I focus on the radio bands which use measurements from the VLA, ATCA, MeerKAT, AMI, LOFAR, ATA, SMA, ASKAP, e-MERLIN, and NEOMA telescopes, which are present in combinations for most raw datasets I use. Most of these have distinct frequency bands of 1.3, 2.4, 4.9, 8.5, and 15 GHz. The frequencies I observe in the project are 1.3, 4.9, 8.5, and 15 GHz due to the higher quality of data taken there. 

## Pipeline
The project relies heavily on the pandas library, mainly for dataframe operations which are critical for calculating the spectral index. Matplotlib and numpy are also used. Custom functions I built to clean, filter, analyze, and plot the data are also used which are defined and explained in data/analysis/grb_functions.py. To operate this pipeline, the filename of the GRB dataset, the redshift of the GRB, and the luminosity distance of the GRB are needed. Before including each dataset in the pipeline, I research the GRB to determine those astronomical values and include a reference paper for future use, which is detailed in the file GRB_info.txt

Each GRB is plotted and briefly analyzed in its separate script in order to evaluate the quality of data before downloading and analyzing spectral data. The scripts are named by {GRB name}.py in src/individual. The pipeline first passes the filename into the pandas read csv function with inferred column names based on known dataset formats, which I edit if need be. I define the redshift and luminosity distance values, and convert the flux and frequency columns into the desired units(milliJy, GHz) for future comparison across datasets. Then, I use my custom cleaning functions that remove invalid and outlier flux values from the dataset. Before filtering or operating on the data, I define lower and upper bounds on the radio epochs detailed above that are adjusted for redshift. This is done manually after observing the datasets and their frequency values, which is necessary due to the slight variations between telescope observation frequencies which can't be known before reading in and observing the dataset. Then, I separate the dataset into four subframes for each distinct band(1.3, 4.9, 8.5, 15 GHz) using a custom function. To gain a visual sense of the data, I then plot the light curves(Flux vs time) of each frequency. 

## Usage
Calculating the spectral index is done through another custom function, of which the process is described in detail in the functions script. I calculate the spectral index between 1.3-4.9, 4.9-8.5, and 8.5-15 GHz. After calculating, I plot these values over time as well. Between the light curve and spectral index plots, I can determine on my own discretion(and consulting my advisor) on which datasets have enough spectral data to keep and analyze. If the data is sufficient, I download the spectral data from the pandas dataframe to a csv file into a seperate folder of spectral data, in data/processed. Not all datasets have sufficient data in each band, and the 4.9-8.5 GHz band ends up having the greatest quantity of spectral data.

With spectral data collected, I use spectral_index_comp.py to analyze the combined data in src/analysis. This plots the spectral datasets of GRBs with quality data mapping directly onto observed times. Additionally, to visually compare the temporal evolution more easily, I plot the spectral data with times shifted to map the index where α = 0(approximately center). These plots reveal the differences in evolution between spectral indices corresponding to characteristic frequencies.

## Future Work
Future work for this project will focus on creating a function to mathematically analyze the difference in spectral index evolution, and measure energetics parameters such as the fraction of energy stored in the magnetic field of the jet(Ɛ_B) and the fraction of energy stored in electrons(Ɛ_e). I also plan to rework this pipeline to omit as much manual adjustment of scripts as possible using a dictionary of parameters corresponding to individual GRBs. Additionally, I am also exploring modeling GRB980425, which was the first GRB with a confirmed supernova counterpart, which reinforced the long GRB progenitor theory. The light curves of this GRB have double peaks in 8.5 and 4.9 GHz but not in 1.3 GHz. Through modeling the singular peaks in each frequency, I will determine if the superposition of these peaks can create an accurate model for the entire light curve. 
