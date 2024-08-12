/***************************************************************************
                          floatops.h  -  description
                             -------------------
    begin                : Thu Apr 22 2004
    copyright            : (C) 2004 by Jose Ignacio Garzon
    email                :

		3D Volume Header functions
 ***************************************************************************/

#ifndef _floatops_hpp
#define _floatops_hpp

#include "libtools/include/Surface.h"
#include "vlvoliter_linear.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



 typedef struct // required by convoluteK_nopad_par() (parallel)
 {
	int chunk; // chunk of slices (depends on the number of threads)
 	int border; // kernel size
 	int i_initial; // Initial i-slice index
 	int i_final; // Final i-slice index
 	int j_initial; // Initial j-slice index
 	int j_final; // Final j-slice index
 	int z_initial; // Initial z-slice index
 	int z_final; // Final z-slice index
 	int stepx; // Volume's step-x
 	int stepy; // Volume's step-y
 	int stepz; // Volume's step-z
 	float *container; // Volume pointer
 	float *container2; // Dummy vol. pointer
 	float *kernel; // Convolution kernel pointer
 	pthread_mutex_t *p_mutex_begin; // mutex pointer
 	pthread_cond_t *p_cond_begin; // condition pointer
 	pthread_mutex_t *p_mutex_end; // mutex pointer
 	pthread_cond_t *p_cond_end; // condition pointer
 	int *p_nended; // pointer to the number of "ended" jobs
 	bool begin; // bool telling whether the job started or not
 } convoluteK_data;


///Namespace of operations of Volumes of Float data
 namespace FOPS
 {
  // BASIC MATHEMATIC OPERATIONS
  /**
  * Add the value to all the voxels of a Volume
  *
  * @param vol: Volume over the operation has place
  * @param value: value to be added
  */
  void add(vlVolume *vol, float value);

  /**
  * Add two volumes in the first of them
  * Alert: Both volumes must have the same size
  *
  * @param vol: First volume. The result will be stored in this volume
  * @param vol2: Second volume
  */
  void add(vlVolume *vol, vlVolume *vol2);

  /**
  * multiplies the value to all the voxels of a Volume
  *
  * @param vol: Volume over the operation has place
  * @param value: value to be multiplied
  */
  void mul(vlVolume *vol, float value);
  /**
  * Return the addition of the values of all voxels (only positives)
  *
  * @param vol: Volume over the operation has place
  */
  float calc_total(vlVolume *vol);
  /**
  * Return the addition of the values of all voxels (only positives)
  *
  * @param vol: Volume over the operation has place
  */
  float calc_total_dens(vlVolume *vol);
  /**
  * Return the addition of the values of all voxels
  *
  * @param vol: Volume over the operation has place
  */float calc_total_posneg(vlVolume *vol);
  /**
  * Return the minimum value of the volume
  *
  * @param vol: Volume over the operation has place
  */
  float min_mass(vlVolume *vol);
  /**
  * Return the maximum value of the volume
  *
  * @param vol: Volume over the operation has place
  */
  float max_mass(vlVolume *vol);
  /**
  * Return the average value of the volume
  *
  * @param vol: Volume over the operation has place
  */
  float calc_average(vlVolume *vol);
  /**
  * Return the sigma of all the values of the volume
  *
  * @param vol: Volume over the operation has place
  */
  float sigma(vlVolume *vol);
  /**
  * Return the normal of all the values of the volume
  *
  * @param vol: Volume over the operation has place
  */
  float normal(vlVolume *vol);
  /**
  * Return the maximum distance between the center of the volume
  * (not the center of masses) and the voxel with positive value farest away
  * (in voxels)
  *
  * @param vol: Volume over the operation has place
  */
  int max_lenght(vlVolume *vol);
  /**
  * Return the center of masses of the volume (weighted by the mass of ecah voxel)
  *
  * @param vol: Volume over the operation has place
  * @param cx,cy,cz: returned position of the center (in voxels)
  * */
  void center_masses(vlVolume *vol, float *cx,float *cy, float *cz);
  /**
  * Return the geometric center of the volume (Depending on the voxels with positive value)
  *
  * @param vol: Volume over the operation has place
  * @param cx,cy,cz: returned position of the center (in voxels)
  * */
  void center_geometric(vlVolume *vol, float *cx,float *cy, float *cz);
  /**
  * Return the middle center of the volume
  *
  * @param vol: Volume over the operation has place
  * @param center: returned position of the center (in voxels)
  * */
  void center_vol(vlVolume *vol, float *center);
  /**
  * Return the Maximum Radius from the center of masses to the positive voxels (in voxels)
  * @param vol: Volume over the operation has place
  */
  float radiusMax(vlVolume *vol);
  /**
   * Mon made (6/2/2009)
   * Computes Average inside a frame (size = frame dimensions = padding)
   */
   float avg_frame(vlVolume *vol, vlDim size);
  /**
   * Mon made (13/9/2010)
   * Computes Average (DOUBLE precision) inside a frame (size = frame dimensions = padding)
   */
   double avg_frameD(vlVolume *vol, vlDim size);

  /**
   * Mon made (6/2/2009)
   * Computes Sigma inside a frame (size = frame dimensions = padding)
   */
   float sig_frame(vlVolume *vol, vlDim size);
  /**
   * Mon made (13/9/2009)
   * Computes Sigma (DOUBLE precision) inside a frame (size = frame dimensions = padding)
   */
   double sig_frameD(vlVolume *vol, vlDim size);




  // Filters in the Gaussian Domain
  /**
  * Gaussian filter over a volume
  *
  * @param fourier: Volume
  * @param res: Resolution of the gaussian filter
  * 				Limits: Minimum: Voxel size
  * 						Maximum: voxel size x map lenght /2
  * @param option: Low-pass (0) or high-pass filter (1)
  */
  void GaussFilter(vlVolume *fourier,float res, int option);
  /**
  * Low Gaussian Filter (Pablo's version)
  * @param fourier: Volume
  * @param res: Resolution of the gaussian filter
  * 				Limits: Minimum: Voxel size
  * 						Maximum: voxel size x map lenght /2
  *
  * Pablo's
  *the constant a for a Gaussian y=exp(-a*x^2) is
  *        a=ln(2)/res^2   for res defined at y=0.5
  *        a=1/res^2       for res defined at y=1/e
  *
  *  here we use the res defined in Fourier space
  *  eg. Fourier transform of the above function: Y=exp(-pi^2*k^2/a)
  *   a=ln(2)*pi^2/res^2   for res defined at Y=0.5
  *   a=pi^2/res^2         for res defined at Y=1/e
  *                          reciprocal of the 1/2 width of a Gaussian in Fourier space
  *
  *   Notes:  the gaussian exp 3/2 factor vanishes due to
  *                          sigma_3D=sigma_1D*sqrt(3)
  *                          sigma1D=(res*1/2)/width
  *
  *           res_0.5=1/sqrt(ln2)*res_1/e = 1.20112240878*res_1/e
  */
  void GaussFilter(vlVolume *vol,float res);

  /** Mon made (13/2/2009)
  * Map filtering in Fourier Space.
  * Needs: two "fftw plans" (direct and inverse FFT), a Fourier map, the resolution,
  * and the fourier normalization factor "fttw_cte".
  * Outputs: the filtered map will be placed in the map pointed by the "fftw plan"
  */
  void GaussFilter(void *ppv, void *ppv2, vlVolume *fourier, float res, float fttw_cte);

  /**
  * Band Gaussian filter
  *
  * @param fourier: Volume
  * @param res1: low-Resolution limit
  * @param res2: high-Resolution limit
  */
  void GaussBandFilter(vlVolume *vol,float res1,float res2);
  /**
  * Butter Filter
  *
  * @param vol: Volume
  * @param res: Resolution limit of the filter
  * 				Limits: Minimum: Voxel size
  * 						Maximum: voxel size x map lenght /2
  * @param n: Filter order. The biggest the more accentuated the change in the resolution limit
  * @param option: Low-pass (0) filter or high-pass (1) filter
  */
  void ButterFilter(vlVolume *vol,float res, int n,int option);
  /**
  * Butter Band Filter
  *
  * @param vol: Volume
  * @param res1: Resolution low limit of the filter
  * @param res2: Resolution highlimit of the filter
  * @param n: Filter order. The biggest the more accentuated the change in the resolution limit
  */
  void ButterBandFilter(vlVolume *vol,float res1, float res2, int n);
  /**
  * Laplacian Filter
  *
  * @param vol: Volume
  */
  void laplacianFilter(vlVolume *vol);




  //REDIMENSION OF VOLUME
  /**
  * Padding a Volume (both sides)
  *
  * @param old: Original Volume
  * @margin: extension of the new margins
  * @return New Volumen
  */
  vlVolume * padVolume(vlVolume *old, vlDim margin,bool self=false);
  /**
  * Padding a Volume (only by one side, final side)
  *
  * @param old: Original Volume
  * @margin: extension of the new margins
  * @return New Volumen
  */
  vlVolume * padVolume(vlVolume *old, vlDim margin,int dummy,bool self=false);
  /**
  * Creates a new Volume from an old Volume with a new voxel size
  * It interpolates the values of the new voxels from the old's one
  *
  * @param old: Old volume
  * @param newUnit: New voxel size for the new Volume
  * @param self: The new Volume will be stored in the same pointer of the old one. Delete old volume-
  * @return:  New Volume
  */
  vlVolume *interpolate_map(vlVolume *old, vlUnit newUnit,bool self=false);
  /**
  * Maximum distance betwwen two voxels with positive values in
  * a volume
  *
  * @param vol: Volume
  */
  float maxWidth(vlVolume *vol);

  /**
  * Change the dimensions of a Volume
  *
  * @param old: Volume
  * @param size: New size of the new volume
  * @param self: The new Volume will be stored in the same pointer of the old one. Delete old volume-
  * @return New Volume
  */
  vlVolume * resize(vlVolume *old, vlDim size,bool self=false);


  /* Re-sizes a map defined by its origin and max corner (real space units). (Mon 3/3/2009)
   * The final map will be defined by its "pmin" and "pmax" corners,
   * and will be placed in the same framework as "old" one.
   *
   *  DO IT WITHOUT INTERPOLATION!!! (WHENEVER...)
   *  */
  vlVolume *resize(vlVolume *old, vlPoint3f pmin, vlPoint3f pmax, bool self);


  /**
   * Project the orig volume in the dest volume
   */
  void projectVolume(vlVolume *orig, vlVolume *dest);



  //VOLUME MODIFICATION
  /**
  * Set a minimum for valid values. All voxels with a value below the limit will
  * bet set to 0
  *
  * @param vol: Volume
  * @limit: threshold to be applied
  */
  void threshold(vlVolume* vol, float limit);


  /**
   *  normalize to have the same density at a given cutoff.
   * bet set to 0
   *
   * @param vol: Volume
   * @limit: threshold to be applied
   * @dens: total density
   */
   void norm_threshold(vlVolume* vol, float limit, float dens);


  /**
  * Set a limit for valid values. All voxels with a value closer to zero than the limit will
  * bet set to the zero
  *
  * @param vol: Volume
  * @limit: threshold to be applied
  */
  void threshold_zero(vlVolume *vol, float limit);

  /**
  * Set a limit for valid values. All voxels with a value below the limit will
  * bet set to the limit
  *
  * @param vol: Volume
  * @limit: threshold to be applied
  */
  void thresholdDown(vlVolume *vol, float limit);


  /**
  * Set a limit for valid values. All voxels with a value over the limit will
  * bet set to the limit
  *
  * @param vol: Volume
  * @limit: threshold to be applied
  */
  void thresholdUp(vlVolume *vol, float limit);

  /**
  * Normalize a volume by a factor
  *
  * @param vol: Volume
  * @param factor: Normalization factor
  */
  void normalize(vlVolume* vol, float factor);

  /**
  * Change the values of the voxels to modify the domain of possible values
  *
  * @param vol: Volume
  * @param max: New maximum limit. The voxels with the maximum current value will be set to this value
  * @param min: New minimum limit. The voxels with the minimum current value will be set to this value
  */
  void changeDensity(vlVolume* vol,float max, float min);

  /**
  * Change the values of the voxels to modify the values distribution
  *
  * @param vol: Volume
  * @param mean: New average mean of the distribution
  * @param sigm: New sigma of the distribution
  */
  void changeDensity2(vlVolume* vol,float mean, float sigm);

  /**
  * Find the voxels that form the "surface" of the volume. These are the voxels
  * in contact with empty voxels. Voxels with a value down a limit are considered emptied
  *
  * @param vol: Volume
  * @param limit: Limit to consider a voxel emptied
  * @param boundaryVoxels: Collection of voxels that compound the surface
  */
  void surface(vlVolume* vol, float limit, std::vector<vlPoint3ui> & boundaryVoxels);

  /**
  * Degrades the "surface" of a volume a number of layers
  *
  * @param vol: Volume
  * @param boundaryVoxels: Collection of voxels that compound the surface. MUST be precompute by surface
  * @param iter: Number of layers to be degraded
  * @param newBoundaryVoxels: Collection of voxels that compound the new surface of the volume after the degradation
  */
  void erode(vlVolume* vol, const std::vector<vlPoint3ui> & boundaryVoxels
      ,const uint32 iter, std::vector<vlPoint3ui> & newBoundaryVoxels);

  /**
  * Expands the surface a number of layers
  *
  * @param  vol: Volume
  * @param boundaryVoxels: Collection of voxels that compound the surface. MUST be precompute by surface
  * @param iter: Number of layers to be dilated
  * @param newBoundaryVoxels: Collection of voxels that compound the new surface of the volume after the dilatation
  * @param newTag: Value inserted in the voxels of the new layers
  */
  void dilate(vlVolume* vol, const std::vector<vlPoint3ui> & boundaryVoxels
      ,const uint32 iter, std::vector<vlPoint3ui> & newBoundaryVoxels, float newTag);




  // VOLUME TRANSFORM
  /**
  * Returns the image of the volume in the Fourier Domain.
  *
  * @param vol: Volume in the real domain
  * @return Volume in the Fourier domain
  */
  vlVolume *Fourier(vlVolume * vol);

  /**
  * Returns the image of the real domain volume in the Fourier Domain
  *
  * @param vol: Volume in the real domain
  * @param fourier: Volume in the Fourier domain
  */
  void Fourier(vlVolume * vol,vlVolume *fourier);

  /**
  * Returns the image of the Fourier domain volume in the real domain
  * ATENTION: Output shows disordered quadrants
  * @param vol: Volume
  * @param odd: 0 if the X size is even, 1 if it is odd
  * @return: Real domain Volume
  */
  vlVolume *iFourier(vlVolume* vol,int odd);

  /**
  * Returns the image of the Fourier domain volume in the real domain
  * ATENTION: Output shows disordered quadrants
  *
  * @param vol: Volume
  * @param odd: 0 if the X size is even, 1 if it is odd
  * @param real: Real domain Volume
  */
  void iFourier(vlVolume*vol,int odd,vlVolume *real);

  /**
  * reduces the dimensions of the volume to keep regions that contains
  * values over a threshold.
  *
  *
  * @param vol: Volume
  * @param threshold: Limit. Values down this limit are considered emptied
  * @param self: The new Volume will be stored in the same pointer of the old one. Delete old volume.
  */
  vlVolume *crop(vlVolume *vol, float threshold=0.0,bool self = false);

  /**
  * reduces the dimensions of the volume to keep regions that contains
  * values over a threshold.
  * In addition, it adds the same specified padding to x,y,z dimensions.
  *
  *
  * @param vol: Volume
  * @param threshold: Limit. Values down this limit are considered emptied
  * @param pad: Number of voxels to pad each dimension (+ and - x,y,z)
  * @param self: The new Volume will be stored in the same pointer of the old one. Delete old volume.
  */
  vlVolume *cropad(vlVolume *vol, float threshold=0.0, int pad = 0, bool self = false);

  /**
  * reduces the dimensions of the volume too keep regions that contains
  * values over a threshold Keeping the center of the volume.
  *
  *
  * @param vol: Volume
  * @param threshold: Limit. Values over the positive threshold and down the negative threshold are considered
  * @param self: The new Volume will be stored in the same pointer of the old one. Delete old volume.  */
  vlVolume *crop_square(vlVolume *vol, float threshold=0.00001, bool self = false);
  /**
  * Auxiliar function to refit the cuadrants of a real domain Volume that has been
  * created by a inverse Fourier transform function
  *
  * @param volume: Volume
  */
  void moveCuadrants(vlVolume **volume);


  // REAL DOMAIN FILTERS
  /**
  * Returns a map filtered by a  Kernel. The new volume will be padded in order to allow the introduction of the kernel
  *
  * @param vlm: Volume
  * @param kernel: kernel buffer to apply (created by compute_kernel_Gaussian)
  * @param size: Size of the kernel
  * @param self: The new Volume will be stored in the same pointer of the old one. Delete old volume.
  */
  vlVolume* convoluteK(vlVolume *vlm,float *kernel, int size, bool self=false);

  /**
   * Returns a map filtered by a  Kernel. The volume will NOT be padded.
   * IMPORTANT: The kernel will not be applied over the voxels in the border.
   *
   * @param vlm: Input Volume. IMPORTANT: It must be wide enough to allow introduction of the kernel
   * @param n: Output filtered volume. IMPORTANT: Must be pre-allocated
   * @param kernel: kernel buffer to apply
   * @param size: Size of the kernel
   * @param self: The new Volume will be stored in the same pointer of the old one. Delete old volume.
   */
  vlVolume* convoluteK_nopad(vlVolume *vlm, vlVolume *n, float *kernel, int size, bool self=false);

#ifdef USE_TBB // Enables Intel's Threading Building Blocks parallel routines

  // Convolutes just one "i" slice by a kernel
  void ConvSlice(int i, convoluteK_data const *data);
  // Convolutes one chunk of slices by a kernel
  void ConvSliceX(int i, convoluteK_data const *data);

  // Mon made (3/6/2012)
  // Convolucion de un kernel sobre un mapa using Intel's TBB "parallel_for"
  // vlm= Mapa sobre el que se aplica el kernel
  // kernel= Kernel a aplicar
  // size= Dimension del Kernel
  // Volume must be big enought to allocate the kernel in the border (pad it before!)
  // vlVolume "n" should be already allocated!
  vlVolume* convoluteK_nopad_tbb(vlVolume *vlm, vlVolume *n, float *kernel, int size, int nthreads = 0, int granularity = 0, bool self = false);

#endif

#ifdef USE_PTHREAD // Enables PThread parallel routines
  /**
   * Thread-routine to convolute some cubic region of a map with a given kernel.
   *
   * @param threadarg: Thread parameters
   */
  void *convoluteK_nopad_thread(void *threadarg);

  // Mon made (27/5/2013)
  // Initialization routine for "convoluteK_nopad_par"
  void convoluteK_nopad_par_init(int nthreads, convoluteK_data **p_threads_data, pthread_t **p_threads, float *kernel, int kern_size);

  /**
   * Returns a map filtered by a  Kernel. The volume will NOT be padded.
   * IMPORTANT: The kernel will not be applied over the voxels in the border.
   *
   * @param vlm: Input Volume. IMPORTANT: It must be wide enough to allow introduction of the kernel
   * @param n: Output filtered volume. IMPORTANT: Must be pre-allocated
   * @param kernel: kernel buffer to apply
   * @param size: Size of the kernel
   * @param nthreads: Number of threads used in parallelization.
   * @param self: The new Volume will be stored in the same pointer of the old one. Delete old volume.
   */
  vlVolume* convoluteK_nopad_par(vlVolume *vlm, vlVolume *n, float *kernel, int size, int nthreads, convoluteK_data *threads_data, bool self=false);
#endif

  /**
  * Computes the correlation of two volumes. Both volumes must have the same size.
  *
  * @param vol,vol2: Volumes
  * @return Correlation of both volumes
  */
  float correlation(vlVolume *vol, vlVolume *vol2);

  /**
    * Simple correlation (Multiply vozel * voxel and compute average value).
    *
    * @param vol,vol2: Volumes
    * @return Correlation of both volumes
    */
    float correlation_simple(vlVolume *vol, vlVolume *vol2);

/**
 * Computes and normalizes the correlation inside a cubic-frame (No padding boders) or only using positive voxels
 *
 * @param vol,vol2: Input voxels
 * @param size: kernel width (Determines padded zone)
 * @param nonzero: Inside a cubic-frame (false) using positive voxels (true)
 * @return Correlation of both volumes normalized by the number of used voxels
 */
float correlation_frame(vlVolume *vol, vlVolume *vol2,int size, bool nozero=true);

/**
* Mon made (3/6/2008)
* Cross.Corr. inside a cubic-frame (nozero=false) or inside the final volume (nozero=true, default)
* If volcalc=true then compute BOTH averages & sigmas, otherwise target (vol) avg & sig should be provided.
* size = kernel width
* vol = Target (Final) volume
* If we have avg1 and sig1 (from: vol), we should provide them!
*/
 float correlation_frame(vlVolume *vol2, vlVolume *vol, vlDim size, bool nozero, bool volcalc, float volavg, float volsig);

/**
* Mon made (13/9/2010) (DOUBLE precision in accumulations)
* Cross.Corr. inside a cubic-frame (nozero=false) or inside the final volume (nozero=true, default)
* If volcalc=true then compute BOTH averages & sigmas, otherwise target (vol) avg & sig should be provided.
* size = kernel width
* vol = Target (Final) volume
* If we have avg1 and sig1 (from: vol), we should provide them!
*/
 float correlation_frame(vlVolume *vol2, vlVolume *vol, vlDim size, bool nozero, bool volcalc, double volavg, double volsig);

/** Mon made (3/7/2008)
 * Local Cross.Corr. inside a cubic-frame (nozero=false) or inside the Final volume (vol)
 * (nozero=true, default)--> Warning, only NON-ZERO voxels will be taken into account!
 * size = kernel width
 * vol = Final volume
*/
float localcorr(vlVolume **p_corrmap, vlVolume *vol2, vlVolume *vol,int size, bool nozero);

/**
 *  Mon made (15/7/2008)
 * Volume Spherical Average
 * size = kernel width (to avoid unnecessary computations)
 *radius --> [1,(size-1)/2]
*/
void spherical_avg(vlVolume *vol_in, vlVolume *vol_out, int size, int radius);

/**
 *
 *Mon made (15/7/2008)
 *
 * Weighted Correlation Map Normalization
 * p_norm --> [0,1]
 */
 void norm_wcorr(vlVolume *vol, float p_norm);



/**
 * Computes and normalizes the weighted correlation inside a cubic-frame (No padding boders) or only using positive voxels
 *
 * @param vol,vol2: Input voxels
 * @param volw: voxel weights
 * @param size: kernel width (Determines padded zone)
 * @param nonzero: Inside a cubic-frame (false) using positive voxels (true)
 * @return Correlation of both volumes normalized by the number of used voxels
 */
float correlation_weight(vlVolume *vol2, vlVolume *vol, vlVolume *volw, int size, bool nozero);

/**
 * Creates a map with the average value around a voxel
 *
 * @param vol_in: Input volume
 * @param vol_out: Output volume. It must be pre-allocated
 * @param size:  Kernel width (Determines padded zone)
 * @param radius: Radius around the voxel to compute the average value [1,(size-1)/2]
 * @return Averaged Correlation
 */
void spherical_avg(vlVolume *vol_in, vlVolume *vol_out, int size, int radius);

/**
 * Freak normalization by Mon. Inverts value relevance, the difference between values is set by p_norm
 *
 * @param vol: input volume. It is modified
 * @param p_norm: Set the difference between voxels. All voxels set to 1 (0). Maximum difference (1)
 */
void norm_wcorr(vlVolume *vol, float p_norm);

  /**
  * Creates a kernel for a Gaussian filter of a given resolution
  * Auxiliary function of convoluteK & others
  *
  * @param kernel: Buffer where the kernel will be stored
  * @param dim_vox: size of the kernel (value returned)
  * @param unit_vox: size of the voxels of the volume the kernel is made for
  * @param resolution: Resolution of the Gaussian filter
  * @sigma_factor: determines the width of the kernel
  */
  void compute_kernel_Gaussian(float **kernel,int *dim_vox,float unit_vox,float resolution,float sigma_factor);
  // + laplacian
  void compute_kernel_Lap_Gaussian(float **kernel,int *dim_vox,float unit_vox,float resolution,float sigma_factor);

  /**
  * Mon made (5/2/2009)
  * Dotproduct inside a cubic-frame (avoiding padding...)
  * size = kernel size (avoids non-valid computations)
  * vol = Target (Reference) volume
  */
  float dotprod_frame(vlVolume *vol2, vlVolume *vol, vlDim size);

  /**
  * Returns a mask of indices (array of integers) for fast cross-correlation computations. Mon made (24/8/2018)
  * IMPORTANT: A negative integer indicates array end!
  * vol = Input volume (map)
  * thr = Threshold to consider "in-mask" voxels (values above "thr" will be indexed)
  */
  int *indices_mask(vlVolume *vol, float thr);

  /**
  * Returns a mask of indices (array of integers) with two maps. Values above or below the respective thresholds can be selected.
  * (Mon made 27/8/2018)
  * IMPORTANT: A negative integer indicates array end!
  * vol   = Input mask volume (map)
  * thr   = Threshold to consider "in-mask" voxels
  * over  = If TRUE, values ABOVE "thr" will be indexed, otherwise the values BELOW "thr" will be indexed.
  * vol2  = Input mask volume 2 (map)
  * thr2  = Threshold 2 to consider "in-mask" voxels
  * over2 = If TRUE, values ABOVE "thr2" will be indexed, otherwise the values BELOW "thr2" will be indexed.
  */
  int *indices_mask(vlVolume *vol, float thr, bool over, vlVolume *vol2, float thr2, bool over2);

  /**
  * Create a masked output map from the input "vol" map and an integers-array mask
  * (Mon made 27/8/2018)
  * IMPORTANT: A negative integer indicates array end!
  * vol = Input volume (map)
  * mask = Mask of indices (array of integers), e.g. for fast cross-correlation. WARNING: A negative integer must indicate the array end!
  */
  vlVolume *apply_mask(vlVolume *vol, int *mask);

  /**
  * Returns the dot product between two maps (faster than cross-correlation). Mon made (24/8/2018)
  * IMPORTANT: A negative integer indicates array end!
  * vol  = Input volume 1
  * vol2 = Input volume 2
  * mask = Mask of indices (array of integers) for fast cross-correlation. WARNING: A negative integer must indicate the array end!
  */
  float dotprod_mask(vlVolume *vol, vlVolume *vol2, int *mask);

  /**
  * Returns the normalized cross-correlation between two maps. Mon made (30/8/2018)
  * IMPORTANT: A negative integer indicates array end!
  * vol  = Input volume 1
  * vol2 = Input volume 2
  * mask = Mask of indices (array of integers) for fast cross-correlation. WARNING: A negative integer must indicate the array end!
  */
  float correlation_mask(vlVolume *vol, vlVolume *vol2, int *mask);


  //WRITE
  /**
  * Write a volume in Situs format
  *
  * @param vol: Volume
  * @param: Name of the file to be written
  */
  void write(vlVolume *vol,char *fileName);
  /**
  * Write a volume in CCP4 format
  *
  * @param vol: Volume
  * @param: Name of the file to be written
  * @param flag: Show/Don't show messages
  */
  void writeCCP4MRC(char *filename,vlVolume *vol,bool changendian=false,bool flag=false);

  /**
  * Write a volume in MRC format
  *
  * @param vol: Volume
  * @param: Name of the file to be written
  * @param flag: Show/Don't show messages
  */
  void writeMRC2014(char *filename,vlVolume *vol,bool changendian=false,bool flag=false);




  /**
  * Write a volume in BRIX format
  *
  * @param vol: Volume
  * @param: Name of the file to be written
  * @param flag: Show/Don't show messages
  */
  void writeBRIX(char *filename,vlVolume *vol,bool flag=false);
  /**
  * Write a volume in MATHEMATICA format
  *
  * @param vol: Volume
  * @param: Name of the file to be written
  * @param indian: indian of the data
  */
  void writeMAT5(char *filename,vlVolume *vol,  char *name,int indian);
  /**
  * Write a volume in the format that the filename extension shows.
  *
  * @param vol: Volume
  * @param: Name of the file to be written. The extensions of this name determines the format of the file
  */
  bool writeFile(vlVolume *vol,char *filename,bool ccp4Endian=false);


  //READ
  /**
  * Read a file in Situs format
  *
  * @param filename: Name of the file to be read
  */
  vlVolume* readVol(char *filename);
  /**
  * Read a file in BRIX format
  *
  * @param filename: Name of the file to be read
  * @fparam flag: Shows/does not show messages
  */
  vlVolume* readBrix(char *filename,bool flag=false);
  /**
  * Read a file in CCP4/MRC format
  *
  * @param filename: Name of the file to be read
  * @param flag: Shows/does not show messages
  */
  vlVolume* readCCP4(char *filename,bool flag=false);
   /**
   * Read a file in MRC 2014 format
   *
   * @param filename: Name of the file to be read
   * @param flag: Shows/does not show messages
   * @param unit: Voxel size of the read volume
   * @param indian: Indian of the data
   */
  vlVolume* readCCP4MRC(char *filename,bool flag=false);
   /**
   * Read a file in MAT5 format
   *
   * @param filename: Name of the file to be read
   * @param flag: Shows/does not show messages
   * @param unit: Voxel size of the read volume
   * @param indian: Indian of the data
   */
  vlVolume* readMAT5(char *filename,float unit,int indian=0);
  /**
   * Read a file in ACNT (Sybyl) format
   *
   * @param filename: Name of the file to be read
   */
  vlVolume* readACNT(char *filename);

  /**
  * Read a file in format specified by the extension of the file
  *
  * @param  filename: Name of the file to be read. The extension must match the format of the file
  * @param unit: Voxel size of the read volume (if necessary)
  */
    vlVolume *readFile(char *filename,float unit=1);


  //ROTATIONAL METHODS
  /**
   * Auxiliary function. Applies a a and b dependent shear  over x-axis x of volume
   */
  bool xshear(vlVolume *vol, float a, float b, float s);
  /**
   * Auxiliary function. Applies a a and b dependent shear  over y-axis x of volume
   */
  bool yshear(vlVolume *vol, float a, float b, float s);
  /**
   * Auxiliary function. Applies a a and b dependent shear  over z-axis x of volume
   */
  bool zshear(vlVolume *vol, float a, float b, float s);
  /**
   *Auxiliary function. Interpolates the values of two columns in a rotation
   */
  bool shifter(int n,int nup,float af,float *f,float ag,float *g);
  /**
   * Rotates a volume. It does not work properly for all the rotations (introduce distorsions)
   *
   * @param vol: Volume to be rotated
   * @param psi,theta,phi: Euler angles
   * @param sh1,sh2,sh3: Movement in xyz (in voxels)
   * @param opcion: Euler angle convention. 0 (ZXZ) 1 (XYZ) 2 (ZYZ)
   */
  void rotate(vlVolume *vol,float psi, float theta, float phi,float sh1,float sh2,float sh3,int opcion=0);
  //Funciones de Flip
  /**
   * Auxiliary function. Rotates 180º over z-axis
   */
  bool flipxy(vlVolume *vol);
  /**
   * Auxiliary function. Rotates 180º over z-axis
   */
  bool flipyz(vlVolume *vol);
  /**
   * Auxiliary function. Rotates 180º over y-axis
   */
  bool flipxz(vlVolume *vol);
  /**
    *Rotation of a volume. Direct interpolation
    *
    *@param vol: Original volume
    *@param       *@param vol2:rotate volume
    *@param
    *@param psi,theta,phi: Euler angles
    *@param conv: Convention of Euler angles. 0 (ZXZ) 1 (XYZ) 2 (ZYZ).
    *@return New rotated volume
    */
void rotate_interp(vlVolume *vol, vlVolume *vol2,float psi, float theta, float phi, int conv);

  /**
      *Rotation of a volume. Direct interpolation
      *
      *@param vol: Original volume
      *@param psi,theta,phi: Euler angles
      *@param conv: Convention of Euler angles. 0 (ZXZ) 1 (XYZ) 2 (ZYZ).
      *@return New rotated volume
      */
vlVolume* rotate_simple(vlVolume *vol, float psi, float theta, float phi, int conv);



  //Masks
  /**
  * Creates a masks from a volume. A mask is a volume with values 0 or 1.
  * The voxels with value 1 are voxels in the original volume with a value over a
  * cutoff, the voxels with 0 have a value below the cutoff in the original volume
  *
  * @param vol: Original Volume
  * @param cutoff: Limit to determine if a voxel is valid or invalid
  * @return Mask volume
  */
  vlVolume* createMask(vlVolume *vol, float cutoff);
  /**
  * Returns true if the given position is a valid voxel in a Volume Mask,
  * false in other case
  *
  * @param mask: Volume Mask
  * @position: Position of the voxel to check
  * @return true if the volume is filled, otherwise false
  */
  bool inMask(vlVolume *mask,const vlPoint3i & position);
  /**
  * Reduces the surface of the mask a number of layers
  *
  * @param mask: Volume Mask
  * @param dist: Number of layer to erase
  */
  void eatMask(vlVolume *mask, uint32 dist);
  /**
  * Increases the surface of the mask a number of layers
  *
  * @param mask: Volume Mask
  * @param dist: Number of layer to increase
  */
  void dilateMask(vlVolume *mask, uint32  dist);
  /**
  * First increases and then reduces the mask a number of layers.
  * It provides a softer mask
  *
  * @param mask: Volume Mask
  * @param dist: Number of layer to increase/Decrease
  */
  void beatMask(vlVolume *mask, uint32 layers);
  /**
  * Invert the voxels of a mask. The valid voxel are invalid and viceversa
  *
  * @para, mask: Volume Mask
  */
  void invert(vlVolume *mask);
  /**
  * Fill holes of size fill ?????
  */
  void fill_mask(vlVolume *mask, int fill);
  /**
  *  Count voxels with value over the cutoff
  */
  int count_points(vlVolume *vol, float cutoff);
  /**
  *  Fills a cube of a given length centered in a given postion
  *
  *  @param mask volume
  *  @param cx,cy,cz: Position to center the cube
  *  @param sampling: radius of the cube
  */
  void fillcenter_onmask(vlVolume *mask, float cx, float cy, float cz, float sampling);
  /**
   * creates a map with regions of voxels with a value over a cutoff. A region is created by voxels connected with values over the cutoff.
   * The voxels of each regions are marked with different values
   * It also returns a list indicating, for each region, the number of voxels.
   * It is also possible to delete the smallest regions
   *    *
   *@param vol: Initial map
   *@param cutoff: connected voxels with values over this cutoff will be considered in a same region
   *@param list: array with the number of voxels in each region created
   *@param num_max: Number of created regions
   *@param min_voxels: Minimum number in regions to be considered wide enough
   *@param erase: if true, deletes a region if it is not considered wide enough
   *@param connectivity: Type of connectivity between voxels. 0: minimum connectivity (only 6 neighnors) 1: medium connectivity (22 neighbors) 2: high connectivity (26 neighbors)
   */
  vlVolume * regions(vlVolume *vol, float cutoff, int **list, int *num_max,int min_voxels=-1, bool erase=false ,int connectivity=0);
  /**
   * Labels in a new map all voxels over a cutoff which conform a region with an initial voxel
   *
   *@param vol: Initial map
   *@param mask: result map. It could be filled with previous regions. If the initial point is inside one of this regions, the function does nothing.
   *@param initial_point: initial voxel for the creation of the region
   *@param mark: value to identify the voxels in the region in the result map.
   *@param cutoff: connected voxels with values over this cutoff will be considered in the region
   *@param min_voxels: Minimum number in the region to be considered wide enough
   *@param erase: if true, deletes the region if it is not considered wide enough. To erase the region is marked with -1
   *@param connectivity: Type of connectivity between voxels. 0: minimum connectivity (only 6 neighnors) 1: medium connectivity (22 neighbors) 2: high connectivity (26 neighbors)
   */
  int flood(vlVolume *vol, vlVolume *mask, vlPoint3ui initial_point, int mark, float cutoff, int min_voxels=-1, bool erase=false, int connectivity=0);



  //SURFACE METHOD
  /**
  * Auxiliar function to compute Surface
  * Project the volume in a Surface Grid Representation
  */
  vlVolume *projectSurface( vlVolume *vol,float PR, bool inner=false);

  /**
  * Auxiliar function to compute Mask without holes
  * Project the volume in a Surface Grid Representation
  */
  //vlVolume *surfaceMask( vlVolume *vol,float PR, float PR2,bool inner=false );
  S_Grid surfaceMask_grid( vlVolume *vol,float PR, bool inner=false );
  vlVolume *erode_grid( vlVolume *vol,S_Grid grid,float PR,int offset);
 }

 namespace IOPS{
   /**
   * Auxiliar function to refit the cuadrants of a real domain Volume that has been
   * created by a inverse Fourier transform function
   *
   * @param volume: Volume
   */
   void moveCuadrants(vlVolume **volume);
 }

 #endif
///////////



