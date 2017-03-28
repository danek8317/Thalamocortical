======== README from NEURON version ====================== 

README from Version of the code
https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=82894 For
changes and modification to this code, **SEE BELOW**

NEURON version of Traub et al J Neurophysiol. 2005 Apr;93(4):1829-30.
Single-column thalamocortical network model exhibiting gamma
oscillations sleep spindles, and epileptogenic bursts.

See:
http://senselab.med.yale.edu/senselab/modeldb/ShowModel.asp?model=45539

Prepare for running with nrnivmod mod

See the comparison of each cell type with the fortran output with
nrngui onecell.hoc

Selecting the "Exact traub style ri" which forces the NEURON
connection coefficients to be exactly the same as computed by the
Traub's equivalent circuit style from the FORTRAN demonstrates that
cell and channel properties are reliably represented in this NEURON
translation of the FORTRAN code.

Reliability:

Since the NEURON and Fortran models do not produce quantitatively
identical results there is always some question as to whether
simulation differences are due to substantive parameter translation
errors or can be attributed to different numerical methods.  It must
be realized that our experience has been that every test into a new
runtime domain has exhibited discrepancies that were ultimately
resolved by fixing translation errors.

And also that our comparison tests are only between NEURON and an
already significantly modified FORTRAN code. The bulk of the FORTRAN
modifications are toward more generic FORTRAN syntax to allow the
original ifc compatible FORTRAN to run under g77. Most of the
execution places where the g77 version differs from the ifc version
are straight forward transformations of bulk array assignment into
equivalent elementwise assignment via do loops.  Did we get them all?
Did we assign over ALL the elements in each array?  We did manually
review all ifc to g77 editing changes but a few cases involved our
judgment with regard to whether there was a bug in the original ifc
fortran version.  The modified FORTRAN used for the NEURON comparisons
is available from this model=45539 page. As is, the g77 FORTRAN model
can only be run as 14 processes, one for each cell type and a full
model run takes 20 hours or so. Simplifying to 1/10 the number of
cells gives a model that takes approximate 1.5 hours for 100 ms of
simulation time. Our last network bug, based on significant spike
raster plot discrepancies, was found using a 10 ms run.

We consider the translation of the 14 individual cell types to be
quite reliable based on the quantitative similarity of the g77 and
NEURON isolated 100 ms cell trajectories at the spike detection
location for 0 and 0.3 nA constant current stimulation into the
soma. Note that quantitative similarity demands compartment coupling
of exactly the same values used by the FORTRAN version algorithms
(imitated in NEURON using the "traub_exact()" algorithm where some
branch points had the form of "wye" equivalent circuits, some had the
"delta" form, and all had a different view of how resistance from
child to parent should be computed.)

Network topology and chemical synapse parameter reliability is limited
to the diagnostic power of our specific tests. For quantitative
comparisons we printed the precise FORTRAN network topology to files
and used that information to define the NEURON network
connections. For 10 ms with a 1/10 size network we focused on
quantitative similarity of the spike raster plots. The FORTRAN version
has a spike resolution time of 0.1 ms and all synaptic conductance
trajectories are step functions with that resolution (the underlying dt
is 50 times smaller, dt = 0.002). We prepared a special version of the
NEURON executable to force spike threshold detection on 0.1 ms
boundaries to allow convenient comparison of spike rasters. For the
first 10 ms we judged whether spike discrepancies were due to
FORTRAN-NEURON spikes straddling the 0.1 ms boundaries or whether the
discrepancy was likely to be due to a topology or synaptic parameter
error. The judgment was based on the details of the voltage
trajectory at the spike detector compartment.  We believe that careful
analysis of the first 10 ms of the 100 ms spike raster overlap plot
for the FORTRAN (fat red marks) and NEURON (thin black marks) in
combination with the spike trajectory sensitivity of suppyrRS cells
with respect to number of spikes in their burst after the first spike
will convince that we have gone as far as possible with quantitative
spike location similarity as a diagnostic technique. Further
diagnostics will likely have to be based on specific questions in
regard to qualitative discrepancies and a focus on the NEURON model
itself as the tool for exploration in terms of certain properties
added or subtracted from successive runs. Unfortunately that can only
be done in response to a specific suspicion on the part of the
user. Enumerated below are the known discrepancies between the
representations of the Traub model in FORTRAN and NEURON:

The NEURON nmda saturation is turned off. See the NMDA_saturation_fact
in the FORTRAN groucho.f file and the nrntraub/mod/traub_nmda.mod
file.  Warning: NEURON will not mimic the FORTRAN merely by setting
the factor to 80.

The groucho.f axon_refrac_time is normally set to 1.5. Our
quantitative tests temporarily set this parameter to 0.5. The NEURON
spike detection algorithm defines a spike as a positive going
transition past the trigger value.

Enumerated below are those major components (of which we are aware)
that are in the model but have not been tested in terms of their
quantitative equivalence to 
FORTRAN: gap junctions.  
long term nmda properties and effects.  
ectopic spikes random current stimulation

The bottom line: The spike rasters for a full g77 FORTRAN run (gap
junctions and current stimulation present) and a full NEURON run with
its own independent random variables (no "traub_exact" connection
coefficients, dt = 0.025 (ms), secondorder = 2, and spike detection
with dt resolution) is presented.

========================== Modifications ==========================

Modifications to the Traub's model in Neuron to record transmembrane
currents Also included is file to generate 3D morphology See
Thalamocortical_imem folder

Traub model + additional module to record transmembrane currents 
To generate 3d information run 
	generate3dinformation.hoc 

from Thalamocortical_imem/nrntraub_imem/mylib/shape3d/
(e.g nrinvmodl mod  
	  x86_64/special mylib/shape3d/generate3dinformation.hoc)

========Specific changes in the code ========

in init.hoc:

   default_var("spike_compress", 0)
   instead                       
   default_var("spike_compress", 5)

just before running the simulation
     if (save_currents){init_all_data()}

in parlib.hoc

in prun
	while (t < tstop){
		if (pc.id == 0){print t}        
		take_currents()
		order_write_data()	
		fadvance()
                	}
instead 

in manage_setup.hoc
   just after {localloadfile("hoc/traubcon.hoc")}:
         {load_file("mylib/lib.hoc")}
         {load_file("mylib/params.hoc")}
        
   just before {localloadfile("net/groucho.hoc")}:
         {load_file("mylib/stimulus/stimulus.hoc")}

   just after {localloadfile("net/groucho.hoc")}:
         {load_file("mylib/stimulus/set_stimulus.hoc")}
         {load_file("mylib/init.hoc")}

in mod folders additional files:

   izap.mod
   izap2.mod
   vecevent.mod

Additional mylib library

* init.hoc 
  * initilizes all needed objects like conteners to save
voltage, transmembrane currents e.c.t.  
  * close some or all of the ion channels if passive_axon = 1,
passive_soma = 1, passive_dendrites = 1, passive passive_axon = 1 or
set_zero_fast_na = 1
 
* lib.hoc 

  * contains set of funcions used in other files

params.hoc 

  * set of prameters used to save data, to apply additional stimulus,
to close some/all of ion chanels

modules: 
* saving_data 
  * saving_data.hoc - contains definition of
functions to get, interpolate (with a constance time step -> "t_step"
parameter in params.hoc - since simulation are with a variable time
step) and write the data
			
* additional_curr.dat 

  * specification of single currents which you want to write, number
of lines depends on "n_currents" and "n_not_ionic_curr" in a
"params.hoc"

   1."n_not_ionic_curr" lines which contains name of not ionic current
with sign {-1,0,1} which tells how the current confluence
transmembrane current

   2."n_currents" lines with conditions when the current should be
   truck, e.g you can't truck Ca current if there are no Ca channels in
   the segment

   3."n_currents" lines with definition of currents which you want to
   track ss objects is defined in saving_data.hoc contains
   not_ionic_sources defined at the top of the file

* shape3d 

  * cells_shape _shape.dat - single cell 3d morphology

  * generate3dinformation.hoc - code in hoc to generate random 3d
  positions

  * extra_shape.hoc -

  * shapes_mat.txt - mapping from the population to morphology type
    from cells_shape folder (e.g. n-th line refer to the n-th
    population population)

* stimulus 
  * stimulus.hoc - definition of procedures to set constant or
sinusoidal stimulus (also stimulating the network by virtual
electrode, but not sure if that works) vcevent.hoc - procedures to
create stimulator which stimulates network by artificial spikes, path to
the file with artificial spike times is defined in "spikes_dir" in
params.hoc The stimulator is created if "use_vecstim = 1" is in
params.hoc set_stimulus.hoc - set stimulus (constant or sinusoidal) -
execute procedures defined in stimulus


========================= Running simulation =========================

1. Modify parameters, set stimulus in files
	 init.hoc, mylib/params.hoc and mylib/set_stimulus.hoc

2. Compile modfiles:
	nrnivmodl mod
	( folder i686, x86_64 , powerpc or other will appear)

3. To run the simulation

   mpirun -np n_cores i686/special -mpi init.hoc 
   (depending on an architecture there can be x86_64/special or powerpc/special e.t.c)

   To calculate 3d positions:
	
   mpirun -np n_cores i686/special -mpi generate3dinformation.hoc

====================== TODO ======================

n_currents = 0 creates file, only doesn't write the currents (but
write time , cell, segment e.t.c)

===================== Convert to NSDF =====================

See folder convert_nsdf

Step 1: sh reformat.sh

IMPORTANT, this assumes, the output files from simulation are in
ascending order of time.  pushes the 5th column element into a file
named by the second column element; better organized data and requires
GNU parallel for this step.

Creates folder, ./i/tmp for currents and ./v/tmp for potentials.
After this ./i/*.dat and ./v/*.dat can be deleted because they are
redundant.

Step 2: python convert2hdf5.py

Looks for the ./i/tmp and ./v/tmp folders and subfolders and
re-processes them into sensible NSDF structures. First pass because
there are some exceptions which are corrected in the next steps, and
there is no meta data (and morphology)

Step 3: python correct_files.py

(accordinly use python correct_files_small.py for 10% models) makes
corrections in the NSDF files for the exceptions. Morphology is added
at this stage

Step 4: python meta_data.py

Adds the meta data of when the simulations were performed, by who,
ownership, licenses, etc.

Step 5: python meta_data_correction.py

Corrects the title names etc for the meta data. These are missed in
the first pass.

After this a sane valid NSDF file is born

================ Analysis scripts ================ 

See folder 'analysis_scripts'

These Python scripts are specific for these simulations, although they can
be extended to generic NSDF files.


1.  lfp_parameters.py

    Has the option of selecting 1D, 2D or 3D electrodes next to the
    cortical column. Assigns default file names for dumping when
    LFP is calculated. Includes cell populations and the model size to
    be used for computing the LFP in the next steps. Also, the dataset to
    be used is assigned here to variable 'h'

    Here is a toggle for selection of 1D or 2D electrode selection for
    subsequent calculations.

2.  calc_lfp.py

    Computes the extracellular potentials using the transmembrane currents
    from the population of cortical cells. The electrode positions and the
    populations to be changed for this computation are taken from
    lfp_parameters.py file. Also included is a function to convert the LFP
    calculated here into a NEO object list, which can be used e.g. in elephant.

    In this script we use the point source approximation, i.e., the compartments are
    treated as point sources placed at the mid point of the
    compartments. It also has the low pass filter function used to compute
    the LFP from extracellular potentials (filtered using 2nd order
    Butterworth filter).

3.  create_plot.py

    Shows the measured potentials in a plot for 1D and 2D electrode
    positions. It also shows the mid points of the compartments used in
    the LFP computation, marks the electrode positions, 
    and shows interpolated potentials recorded using the 1D and 2D
    probes. The plot displayed for 1D case shows potentials versus time (x-axis), like
    for a laminar probe, and for the 2D it is like a plane of electrodes as in MEA at time
    point of 110.5 ms - this can be changed here by the user.

4.  raster_spikes.py

    Shows the spiking activity of the network color coded. 
    The different colors represent different neuron types. Up
    arrow indicates an excitatory neuron. Down represents an inhibitory neuron.

================= figures =====================

See folder 'figures'

These Python scripts use the data provided and generate the figures in
the manuscript.

1. figure1.py 

   Places electrodes in 2D and on the left shows the positions of the 
   mid-points of the compartments in black plus the electrode positions. 
   An inset shows the potentials recorded at the marked electrodes.

   On the right interpolated potential from the recorded electrodes 
   is shown at the time point indicated in the inset figure

2. figure2_twocell_data.py and figure2_twocell_map.py

   Shows the internal structure of the NSDF file from these simulations.

3. figure3.py

   Displays some results of analysis of a small network, including
   A) raster plot, B) LFP, and contributions to LFP (which is extracellular 
   potential filtered below 100 Hz) from individual currents, which come 
   specifically from C) NMDA & AMPA (recorded together) D) GABA
   E) capacitive, F) potassium, G) passive, H) calcium, I) sodium J) two
   kinds of calcium low threshold T type currents not causing [Ca2+]
   influx K) anomalous rectifier, L) all other currents, such as ectopic
   currents and depolarizing currents.

4. figure4.py

   Displays contributions to the extracellular potential from
   different types of current sources. Possible thanks to the different
   datasets presented here.



