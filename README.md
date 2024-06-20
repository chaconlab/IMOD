# IMOD

iMOD is a versatile toolkit to perform Normal Mode Analysis (NMA) in internal coordinates (IC) on both protein and nucleic acid atomic structures. Vibrational analysis, motion animations, morphing trajectories, and Monte-Carlo simulations can be easily carried out at different scales of resolution using this toolkit. You can also have access to the latest version from our <a href="http://imods.chaconlab.org/">iMODS</a> online server

#### References

- López-Blanco JR, Garzón JI, Chacón P. (2011). iMod: multipurpose normal mode analysis in internal coordinates. Bioinformatics. 27 (20): 2843-2850.<a href="http://www.ncbi.nlm.nih.gov/pubmed/21873636"><img src="https://chaconlab.org/images/publications/pubmed.jpg" alt="" align="top" border="0" /></a><a href="https://chaconlab.org/PDF/Bioinformatics2011.pdf"><img src="https://chaconlab.org/images/publications/acrobaticon4.gif" alt="" border="0" /></a>

## User guide
Here we give a brief overview of the necessary commands to use iMOD, but we strongly encourage to follow the tutorials. Right now, there are three different executables:
<ul>
<li><a href="#iMODE">iMODE</a> to obtain IC Normal Modes.</li>
<li><a href="#iMOVE">iMOVE</a> to animate IC Normal Modes.</li>
<li><a href="#iMODVIEW">iMODVIEW</a> to visualize Normal Modes.</li>
<li><a href="#iMC">iMC</a> to perform a Monte-Carlo simulation.</li>
<li><a href="#iMORPH">iMORPH</a> to perform Morphing.</li>
</ul>
This user guide describes the usage of these iMOD components.

### iMODE - obtaining IC Normal Modes

To obtain the IC modes, enter the following command at the prompt.

```
> imode <pdb>;
```
where:
<table class="text" style="width: 550px;" border="0" cellspacing="4" cellpadding="2">
<tbody>
<tr>
<td bgcolor="#f0f0f0" width="85">pdb</td>
<td>PDB input file (required)</td>
</tr>
</tbody>
</table>
<p>The default output is:</p>
<ul>
<li><i>imode.log </i>--&gt; Log-file.</li>
<li><i>imode_model.pdb </i>--&gt; Used PDB model.</li>
<li><i><b>imode_ic.evec </b></i>--&gt; IC Normal modes file.</li>
</ul>
<p>To enable deformability computations use the −d option. In this case, additional files will be saved:</p>
<ul>
<li><i>imode_def.pdb </i>--&gt; PDB file with deformability data in B-factor column.</li>
<li><i>imode_mob.pdb </i>--&gt; PDB file with mobility data in B-factor column.</li>
<li><i>imode_defmob.txt </i>--&gt; Plain text file with deformabiliy, mobility and B-factor data.</li>
</ul>

#### Basic Options

Just add the following basic options to customize the IC modes generation after the minimum command line shown above.</p>
<table class="text" style="width: 550px;" border="2" cellspacing="4" cellpadding="2">
<tbody>
<tr>
<td bgcolor="#f0f0f0" width="80">−h</td>
<td>Displays usage information and exits.</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−m &lt;int&gt;</td>
<td>Coarse-Grained model: 0=CA, 1=C5, 2=Heavy-Atom (default=2).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−o &lt;string&gt;</td>
<td>Output files basename (default=imode).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">-d</td>
<td>Turn on deformability calculations (default=disabled).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">-r &lt;float&gt;</td>
<td>Randomly fixed ratio of Dihedral Coordinates (default=disabled).<br /> Example: 0.7 = 70% of dihedrals will be randomly removed.<br /> Rotational/translational coordinates always mobile.</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">-f &lt;string&gt;</td>
<td>ASCII file defining the ICs to be fixed with the format:<br />
<p>Protein: "n phi chi psi"<br /> NAcid: "n alpha beta gamma chi epsilon zeta"<br /> Inter-chain: "n 6D"</p>
<p>Where "n" is the residue index (0,1,..) and the coordinate name (phi, psi, etc...) can be set to 0(fixed) or 1(mobile). Each one of the 6 inter-chain variables should be specified on separate lines in the following order: x,y,z,Rx,Ry,Rz. Note "n" is just the sequential residue index (starting with 0) and NOT the PDB's residue index.</p>
<p>A demo file can be generated using the --save_fixfile option.</p>
</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">-P &lt;int&gt;</td>
<td>Pairwise interaction potential: (default=0)<br />
<p>0= Sigmoid function (= k/(1+(x/x0)^p), if x &lt; c, else k=0).<br /> 1= Tirion's cutoff (= k, if x &lt; c, else k=0).<br /> 2= Hinsen's function.<br /> 3= Topology &amp; Secondary Structure (--func is mandatory).<br /> 4= edNMA formalism (CA-model only).</p>
<p>By default an extra torsional potential will be added.</p>
</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">-K &lt;string&gt;</td>
<td>Force constants ASCII file with 3 cols.: <br /> Where <i> are the corresponding atomic indices (1,2,...)<br /> A demo file can be generated using --save_Kfile option.</i></td>
</tr>
<tr>
<td bgcolor="#f0f0f0">-n &lt;int/float&gt;</td>
<td>Used modes range, either number [1,N] &lt;integer&gt;, or ratio [0,1) &lt;float&gt; (default=20).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">-x</td>
<td>Considers first CHI dihedral angle (default=disabled).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">-S &lt;string&gt;</td>
<td>All dihedral coordinates with a given secondary structure (SS) will be removed (see --ss).<br /> Example: "HE" will fix the dihedrals corresponding to alpha-helices and beta-sheets.</td>
</tr>
</tbody>
</table>

#### Advanced options

<p>Only for real expert users!</p>
<pre>   --ss &lt;string&gt;
      Secondary Structure ASCII file with 2 cols.: &lt;n&gt; &lt;char&gt;
     Where &lt;n&gt; is the corresponding residue index (0,1,...), and &lt;char&gt; is
     the single character SS identifier. By default SS will be computed
     internally (H=helix, E=strand, C=coil). 
   --save_fixfile
      Save fixation file as &lt;basename.fix&gt; (to be used with -r or -S
     options; otherwise a fully mobile file will be generated)
     (default=disabled). 
   --save_cart
      Save Cartesian modes as &lt;basename_cart.evec&gt; (default=disabled) 
   --save_wcart
      Save Mass-weighted Cartesian modes as &lt;basename_wcart.evec&gt;
     (default=disabled). 
   --save_Kfile
      Save atom-pairwise force constants file as &lt;basename_Kfile.dat&gt; (to
     be used with -K option) (default=disabled).
   --save_SSfile
      Save secondary structure file as &lt;basename.ss&gt; (to be used with -S or
     -P=2 options) (default=disabled). 
   --save_covar
      Saves the predicted covariance matrix at selected Temperature in
     binary packed storage format as &lt;basename_covar.evec&gt;. If --save_wcart
     selected, then mass-weighted covariance matrix will be computed
     instead (default=disabled). 
   --k0_c &lt;float&gt;
      Sigmoid function distance cutoff (default=10A). 
   --k0_k &lt;float&gt;
      Sigmoid function stiffness constant (default=1.0). 
   --k0_x0 &lt;float&gt;
      Sigmoid function inflexion point (default=3.8A). 
   --k0_p &lt;float&gt;
      Sigmoid function power term (default=6). 
   --k1_c &lt;float&gt;
      Tirion's method distance cutoff (default=10A). 
   --k1_k &lt;float&gt;
      Tirion's method stiffness constant (default=1.0). 
   --k2_c &lt;float&gt;
      Non-bonding distance cutoff applied to --func option (default=10A). 
   --nomodel
      Disables PDB model building. Warning: introduced PDB model must match
     the CG selected with the -m option (default=disabled). 
   --nomass
      Disables mass weighting (default=disabled). 
   --notors
      Disables extra torsional potential (default=disabled). 
   --norm
      Enables (norm=1) eigenvector normalization. Note this does not change
     vector direction (default=disabled). 
   --func &lt;string&gt;
      ASCII file defining the force constant functions to be applied
     according to Topology and/or Secondary Structure. The 5 cols. format
     is: &lt;SS&gt; &lt;t&gt; &lt;k&gt; &lt;x0&gt; &lt;pow&gt;
     Where &lt;SS&gt; is the two character pairwise interaction identifier, &lt;t&gt;
     is the topology, and &lt;k&gt;,&lt;x0&gt;,&lt;pow&gt; are the corresponding sigmoid
     function parameters (see -P option). If --ss is not specified, the XX 
     pairwise interaction identifier must be introduced. This way, only 
     topologies will be considered. If &lt;t&gt; is "-1", any previously 
     not-matched topology will be considered.
   --model_out &lt;int&gt;
      Output Coarse-Graining model: 0=CA, 1=C5, 2=Heavy-Atom
     (default=disabled). 
   --chi_out
      Considers first CHI dihedral angle in output modes
     (default=disabled). 
   --save_covar_out
      Computes and Saves the predicted covariance matrix for the output
     model at selected Temperature in binary packed storage format as
     &lt;basename_covarf.evec&gt;. If --save_wcart selected, then mass-weighted
     covariance matrix will be computed instead (default=disabled). 
   -T &lt;double&gt;,  --temperature &lt;double&gt;
      Temperature [K] for covariance matrix computation (default=300). 
   --seed &lt;unsigned&gt;
      Pre-define the random number generator SEED (Mersenne Twister)
     (default=random-seed from /dev/urandom) 
   --verb &lt;int&gt;
      Verbose level (0=low, 1=medium, 2=high) (default=0).
   --,  --ignore_rest
      Ignores the rest of the labeled arguments following this flag. 
   -v,  --version
      Displays version information and exits. 
   -h,  --help
      Displays usage information and exits. <br /><br /></pre>
      
## iMOVE - animating IC normal modes

To show the vibrational motion of a given mode, just enter the following command at the prompt:
```
>imove <in_pdb> <ptraj> <out_pdb> <int>
```
where:
<table class="text" border="2" cellspacing="4" cellpadding="2">
<tbody>
<tr>
<td bgcolor="#f0f0f0" width="85">&lt;in_pdb&gt;</td>
<td>PDB input file. (required)</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">&lt;ptraj&gt;</td>
<td>Normal Modes input file name (.evec). (required)</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">&lt;out_pdb&gt;</td>
<td>Output Multi-PDB file. (required)</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">&lt;int&gt;</td>
<td>Mode number to be moved (1,2,...,size). (required)</td>
</tr>
</tbody>
</table>
<p>The unique output is the Multi-PDB file named &lt;out_pdb&gt;</p>

#### Basic Options

In this section, the basic options to customize animations are detailed.
<table class="text" style="width: 550px;" border="2" cellspacing="4" cellpadding="2">
<tbody>
<tr>
<td bgcolor="#f0f0f0" width="80">−h</td>
<td>Displays usage information and exits.</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−m &lt;int&gt;</td>
<td>Coarse-Graining model: 0=CA, 1=C5, 2=Heavy-Atom, 3=NCAC, 4=CA-only (default=2).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−x</td>
<td>Considers first CHI dihedral angle (default=disabled).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−c &lt;int&gt;</td>
<td>Number of conformations generated (default=11). It should be an odd number!</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−a &lt;float&gt;</td>
<td>Amplitude linear factor to scale motion (default=2).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−T &lt;float&gt;</td>
<td>Temperature [K] (default=300).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−f &lt;string&gt;</td>
<td>Input ASCII file defining the ICs that were fixed during NMA (default=disabled). If modes were computed removing arbitrary ICs, the user must introduce here the file generated by iMode's --save_fixfile option.</td>
</tr>
</tbody>
</table>

#### Advanced options

Only for real expert users!
<pre>   --cart
      Mandatory if Cartesian modes are supplied, otherwise moving in
     Internal Coordinates (default=disabled). 
   --model_out &lt;int&gt;
      Output Coarse-Graining model: 0=CA, 1=C5, 2=Heavy-Atom. The default
     output model will be that selected with the -m option. 
   --linear
      Enables linear motion instead of sinusoidal bounce
     (default=disabled). 
   --mov &lt;int&gt;
      Motion Type (default=2): 0=K-matrix, 1=V/W-arrays, 2=Simple-Rotations
     , 3=Linear (if Cartesian modes). 
   --chi_out
      Considers first CHI dihedral angle in output models
     (default=disabled) (DEVELOPER's). 
   --fixIC &lt;string&gt;
      Plain-text file defining the fixed Internal Coordinates. Each line
     will contain the index (0,1,...) of the ICs to be removed
     (DEVELOPER's). 
   --,  --ignore_rest
      Ignores the rest of the labeled arguments following this flag. 
   -v,  --version
      Displays version information and exits. 
</pre>

## iMODVIEW - normal modes visualization

An alternative way to visualize a normal mode motion is the arrow representation. To this end type:
```
>imodview <pdb> <ptraj/Kfile> <filename>
```
<p>where:</p>
<table class="text" border="2" cellspacing="4" cellpadding="2">
<tbody>
<tr>
<td bgcolor="#f0f0f0" width="85">&lt;pdb&gt;</td>
<td>PDB input file. Warning: This model must match exactly the one used in modes computation. (required)</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">&lt;ptraj/Kfile&gt;</td>
<td>Input eigenvectors file (ptraj). Warning: Only Cartesian modes allowed.</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">&lt;filename&gt;</td>
<td>Output VMD file. (required)</td>
</tr>
</tbody>
</table>
<p>The unique output file (&lt;filename&gt;) can be loaded into VMD using the following command in VMD's terminal:</p>
<pre>source &lt;filename&gt;
</pre>

#### Basic Options

In this section, the basic options to customize mode visualization are detailed.
<table class="text" style="width: 550px;" border="2" cellspacing="4" cellpadding="2">
<tbody>
<tr>
<td bgcolor="#f0f0f0" width="80">−h</td>
<td>Displays usage information and exits.</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−n &lt;int&gt;</td>
<td>Normal Mode index (1,2,...,size) (default=1)</td>
</tr>
</tbody>
</table>

#### Advanced options

<p>Only for real expert users!</p>
<pre>   --color &lt;string&gt;
      Set color. All VMD colors available (default=white). 
   --op &lt;int&gt;&gt;
      Sets the operation method, 1=arrows, 2=springs (default=1). 
   --max &lt;float&gt;&gt;
      Maximum arrow length [A] (default=10). 
   --thick &lt;float&gt;
      Arrow/spring thickness factor (default=0.05). 
   --level &lt;int&gt;
      Sets the averaging level to compute arrows, 0=atoms, 1=residues,
     2=segments, 3=chains (default=0). 
   --pthr &lt;float&gt;
      Minimum percentual amplitude (from maximum) to show arrows
     (default=0, all arrows). 
   --kthr &lt;float&gt;
      Only those springs with force constants above this threshold will be
     shown (default=disabled). 
   --kthr2 &lt;float&gt;
      Only those springs with force constants below this threshold will be
     shown (default=disabled). 
   --,  --ignore_rest
      Ignores the rest of the labeled arguments following this flag. 
   -v,  --version
      Displays version information and exits. <br />
</pre>

## iMC - performing Monte-Carlo simulations

To carry out a basic Monte-Carlo simulation, enter the following command at the prompt:
```
>imc <pdb>; <ptraj>
```
where:
<table class="text" border="2" cellspacing="4" cellpadding="2">
<tbody>
<tr>
<td bgcolor="#f0f0f0" width="85">&lt;pdb&gt;</td>
<td>PDB input file. (required)</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">&lt;ptraj&gt;</td>
<td>Normal modes input file (.evec), either from NMA or PCA. (required)</td>
</tr>
</tbody>
</table>
<p>The default output trajectory will be named <b>imc.pdb</b></p>

#### Basic Options

In this section, the basic options to customize trajectories are detailed.
<table class="text" style="width: 550px;" border="2" cellspacing="4" cellpadding="2">
<tbody>
<tr>
<td bgcolor="#f0f0f0" width="80">−h</td>
<td>Displays usage information and exits</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−m &lt;int&gt;</td>
<td>Input Coarse-Graining model: 0=CA, 1=C5, 2=Heavy-Atom (default=2).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−o &lt;string&gt;</td>
<td>Output files basename. (default=imc)</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−f &lt;string&gt;</td>
<td>Input ASCII file defining the ICs that were fixed during NMA (default=disabled). If modes were computed removing arbitrary ICs, the user must introduce here the file generated by iMode's --save_fixfile option.</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−p &lt;int&gt;</td>
<td>Finds the optimal energy/stiffness scaling factor to obtain the desired average RMSD (Å) from the initial model (default=disabled).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−−Rg &lt;float&gt;</td>
<td>Filter models by target radius of gyration (default=disabled).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−n &lt;int&gt;</td>
<td>Number of eigenvectors to be employed, either number [1,N] &lt;integer&gt;, or ratio from maximum available [0,1) &lt;float&gt; (default=5).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−c &lt;int&gt;</td>
<td>Number of output conformations (default= 100).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−E &lt;float&gt;</td>
<td>Energy/Stiffness scaling factor (mode energy will be multiplied by this value) (default=1.0).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−i &lt;int&gt;</td>
<td>Number of MC iterations per output structure (default=1000).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−x</td>
<td>Considers first CHI dihedral angle (default=disabled).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−T &lt;float&gt;</td>
<td>Temperature [K] (default=300).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−a &lt;float&gt;</td>
<td>Amplitude linear factor to scale motion (default=1).</td>
</tr>
</tbody>
</table>

#### Advanced options

<pre>   --thr &lt;float&gt;
      Enable filtering by absolute tolerance (default=disabled). 
   --Rmsd &lt;float&gt;
      Filter models by target RMSD (default=disabled). 
   --otraj &lt;int&gt;
      Output trajectory format: 0-Normal-Mode, 1-Multi-PDB, 2-AMBER
     (default=1). 
   --cart
      Mandatory if Cartesian modes are supplied, otherwise moving in
     Internal Coordinates (default=disabled). 
   --seed &lt;unsigned int&gt;
      Set the random number generator SEED (Mersenne Twister)
     (default=random-seed from /dev/urandom). 
   --include_first
      Includes input model as first frame in the Multi-PDB trajectory
     (default=disabled). 
   --var
      Input eigenvalues will be considered as variance (pca), otherwise the
     will be force constants (nma) (default=false). 
   --unweight
      Un-mass-weights the input vectors (default=false). The Mass-weighted
     modes (wcart) will be converted into Cartesian coordiantes (cart)
     (DEVELOPER's). 
   --optf &lt;float&gt;
      Factor to scale the optimal energy/stiffness factor
     (default=disabled) (DEVELOPER's). 
   --rfact &lt;float&gt;
      Agressivity factor (default=7.778) (DEVELOPER's). 
   --mov &lt;int&gt;
      Motion Type (default=2): 0=K-matrix, 1=V/W-arrays, 2=Simple-Rotations
     , 3=Linear (if Cartesian modes) (DEVELOPER's). 
   --fixIC &lt;string&gt;
      Plain-text file defining the fixed Internal Coordinates. Each line
     will contain the index (0,1,...) of the ICs to be removed
     (DEVELOPER's). 
   --verb
      Enables verbose. 
   --,  --ignore_rest
      Ignores the rest of the labeled arguments following this flag. 
   -v,  --version
      Displays version information and exits. <br /><br /></pre>


## iMORPH - performing Morphing

To generate a plausible continuous trajectory between two given conformations enter the following command at the prompt.
```
imorph <initial_pdb>; <final_pdb>
```
where:
<table class="text" style="width: 550px;" border="0" cellspacing="4" cellpadding="2">
<tbody>
<tr>
<td bgcolor="#f0f0f0" width="85">initial_pdb</td>
<td>Initial PDB file. (required)</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">target_pdb</td>
<td>Target PDB file. (required)</td>
</tr>
</tbody>
</table>
<p>The trajectory movie will be automatically named <b>imorph_movie.pdb</b></p>

#### Basic Options

In this section, the basic options to customize your morphing are detailed. Just add them after the minimal command shown above.
<table class="text" style="width: 550px;" border="2" cellspacing="4" cellpadding="2">
<tbody>
<tr>
<td bgcolor="#f0f0f0" width="80">−h</td>
<td>Displays usage information and exits</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−m &lt;int&gt;</td>
<td>Coarse-Grained model: 0=CA, 1=C5, 2=Heavy-Atom (default=2).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−o &lt;string&gt;</td>
<td>Output files basename (default=imorph).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−F</td>
<td>Enables full-atom output models (default=disabled).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−r &lt;float&gt;</td>
<td>Randomly fixed ratio of Dihedral Coordinates (default=disabled).<br /> Example: 0.7 = 70% of dihedrals will be randomly removed.<br /> Rotational/translational coordinates always mobile.</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">-f &lt;string&gt;</td>
<td>ASCII file defining the ICs to be fixed with the format:<br />
<p>Protein: "n phi chi psi"<br /> NAcid: "n alpha beta gamma chi epsilon zeta"<br /> Inter-chain: "n 6D"</p>
<p>Where "n" is the residue index (0,1,..) and the coordinate name (phi, psi, etc...) can be set to 0(fixed) or 1(mobile). Each one of the 6 inter-chain variables should be specified on separate lines in the following order: x,y,z,Rx,Ry,Rz. Note "n" is just the sequential residue index (starting with 0) and NOT the PDB's residue index.</p>
<p>A demo file can be generated using the --save_fixfile option.</p>
</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">-P &lt;int&gt;</td>
<td>Pairwise interaction potential: (default=0)<br />
<p>0= Sigmoid function (= k/(1+(x/x0)^p), if x &lt; c, else k=0).<br /> 1= Tirion's cutoff (= k, if x &lt; c, else k=0).<br /> 2= Hinsen's function.<br /> 3= Topology &amp; Secondary Structure (--func is mandatory).<br /> 4= edNMA formalism (CA-model only).</p>
<p>By default an extra torsional potential will be added.</p>
</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−i &lt;int&gt;</td>
<td>Maximum number of iterations (default=100000).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−n &lt;int/float&gt;</td>
<td>Used modes range, either number [1,N] &lt;integer&gt;, or ratio [0,1) &lt;float&gt; (default=0.1). In any case, the value of --addnevs option will be added.</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−x</td>
<td>Considers first CHI dihedral angle (default=disabled).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−s &lt;float&gt;</td>
<td>Initial amplitude applied to the merge displacement (default=5).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−e &lt;int/float&gt;</td>
<td>Excited modes range, either number [1,nevs] &lt;integer&gt;, or ratio [0,1) &lt;float&gt; (default=0.1).</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−R &lt;float&gt;</td>
<td>Randomly fixed ratio of Internal Coordinates (default=disabled).<br /> Example: 0.7 = 70% of IC will be randomly fixed.</td>
</tr>
<tr>
<td bgcolor="#f0f0f0">−S &lt;string&gt;</td>
<td>All dihedral coordinates with a given secondary structure (SS) will be removed (see --ss). Ex: "HE" will fix the dihedrals corresponding to alpha-helices and beta-sheets.</td>
</tr>
</tbody>
</table>

#### Advanced options

<p>Only for real expert users!</p>
<pre>   --ss &lt;string&gt;
      Secondary Structure ASCII file with 2 cols.: &lt;n&gt; &lt;char&gt;
     Where &lt;n&gt; is the corresponding residue index (0,1,...), and &lt;char&gt; is
     the single character SS identifier. By default SS will be computed
     internally (H=helix, E=strand, C=coil). 
   --conv &lt;double&gt;
      Convergence RMSD threshold (default=0.01Å). 
   --conv_win &lt;int&gt;
      Window length (number of iterations) to estimate convergence
     (default=1000). 
   --delta_save &lt;float&gt;
      RMSD increment to save a new trajectory frame (default=0.5Å). If a
     negative integer value is introduced, a new frame will be saved each
     --delta_save iterations. 
   --k0_c &lt;float&gt;
      Sigmoid function distance cutoff (default=10Å). 
   --k0_k &lt;float&gt;
      Sigmoid function stiffness constant (default=1.0). 
   --k0_x0 &lt;float&gt;
      Sigmoid function inflexion point (default=3.8Å). 
   --k0_p &lt;float&gt;
      Sigmoid function power term (default=6). 
   --k1_c &lt;float&gt;
      Tirion's method distance cutoff (default=10Å). 
   --k1_k &lt;float&gt;
      Tirion's method stiffness constant (default=1.0). 
   --k2_c &lt;float&gt;
      Non-bonding distance cutoff applied to --func option (default=10Å). 
   --func &lt;string&gt;
      ASCII file defining the force constant functions to be applied
     according to Topology and/or Secondary Structure. The 5 cols. format
     is: &lt;SS&gt; &lt;t&gt; &lt;k&gt; &lt;x0&gt; &lt;pow&gt;
     Where &lt;SS&gt; is the two character pairwise interaction identifier, &lt;t&gt;
     is the topology, and &lt;k&gt;,&lt;x0&gt;,&lt;pow&gt; are the corresponding sigmoid
     function parameters. If --ss is not specified, the XX pairwise
     interaction identifier must be introduced. This way, only topologies
     will be considered. If &lt;t&gt; is "-1", any previously not-matched
     topology will be considered. 
   --rediag &lt;float&gt;
      RMSD to trigger NMA (default=0.1Å). 
   --morepdbs
      Saves initial (basename_model.pdb) and final (basename_fitted.pdb)
     models. If the -F option is enabled, the fitted CG-model
     (basename_fitCG.pdb) will be saved (default=disabled). 
   --morermsds
      Enables C-alpha RMSD computation (default=disabled). 
   --notraj
      Disables Multi-PDB trajectory movie output (default=disabled). 
   --nomodel
      Disables PDB model building. Warning: introduced PDB model must match
     the CG selected with the -m option (default=disabled). 
   --nomass
      Disables mass weighting (default=disabled). 
   --notors
      Disables extra torsional potential (default=disabled). 
   --nowrmsd
      Disables Gaussian weighted RMSD (default=disabled). 
   --prob &lt;string&gt;
      Normal mode selection probabitity. Each one of the --nex modes will
     be selected and merged from the --nevs subset according to the
     following probabilities (p): (default=var)
     	plain: constant probability
     	var: proportional to the i-th mode variance, p(i)=
     1/eigenvalue(i)
     	line: lineally decreasing probability, p(i)= 1-i/nevs 
   --addnevs &lt;int&gt;
      Increases --nevs value by --addnevs (default=10). 
   --step2 &lt;float&gt;
      Final amplitude applied to excited modes (amplitude will decrease
     linearly from --step to --step2) (default=step/10). 
   --rand_step
      Random amplitude selection between --step and --step2
     (default=disabled). 
   --nex2 &lt;int&gt;
      Final number of excited eigenvectors, either number [1,nevs]
     &lt;integer&gt;, or ratio [0,1) &lt;float&gt;. Number of excited modes will change
     linearly from --nevs to --nevs2 (default=disabled). 
   --rand_excited
      Random number of excited modes between --nex and --nex2
     (default=disabled). 
   --noscv
      Disable the Scaling Collective Variable method (default=disabled). 
   --norand_weight
      Disables random excited modes weighting (default=disabled). 
   --wrmsd &lt;float&gt;
      Sets the RMSD weighting factor for model alignment. The range is
     [0:1], 0=maximum-weighting 1=no-weighting (default=0). 
   --seed &lt;unsigned int&gt;
      Set the random number generator seed (Mersenne Twister)
     (default=random-seed from /dev/urandom). 
   --time
      Enable clocks (default=disabled). 
   --verb &lt;int&gt;
      Verbose level (0=low, 1=medium, 2=high) (default=0). 
   --,  --ignore_rest
      Ignores the rest of the labeled arguments following this flag. 
   -v,  --version
      Displays version information and exits. 
</pre>
