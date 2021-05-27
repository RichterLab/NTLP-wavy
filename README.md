# NTLP-wavy

To install this project on an ND CRC server, navigate to the desired location
in a terminal and run the following commands:  

---
git clone https://www.github.com/RichterLab/NTLP-wavy  
cd NTLP-wavy/build  
module load intel  
module load mvapich2  
module load netcdf  
make  
cd ../  

---

You can use the copyRun script to create new runs based on the files and
parameters of an old run.  The script automatically updates locations in
params.in and the job name in les.run.  The command syntax is ``` ./copyRun
oldRun newRun -param1 value -param2 value -paramN value ```.  For example, to 
create a new run based on ExRun_Mono, but change the number of particles and 
the number of time steps:  

---
./copyRun ExRun_Mono run1 -tnumpart 9e6 -itmax 100000  

---

Finally, navigate to the newly created run1 directory and queue the run with:  

---
cd run1  
qsub les.run  

---

The les.run script will create a matching set of directories in your scratch365
drive.  This is where you will find the simulation output files. 

Also note the helper scripts ```report``` and ```clean```. The report script 
simply tails the run log to the screen.  The clean script deletes files created
during the run, including all files in the run's scratch directory.
