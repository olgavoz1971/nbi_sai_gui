import numpy as np
from astropy.modeling import models, fitting
from scipy import ndimage

center_method_list = ['barycenter', 'maximum_count', 'maximum_count_filt']
# ------------------ FIT MODELS -------------------------
def fitGauss1D(x,y,a0,m0,stddev0):
 g_init = models.Gaussian1D(amplitude=a0, mean=m0, stddev=stddev0)
 fit_g = fitting.LevMarLSQFitter()
 g = fit_g(g_init, x, y)
 fwhm = 2 * g.stddev.value * np.sqrt(2 * np.log(2))
 return fwhm,g

def fitGauss1D_Const(x,y,a0,m0,stddev0):
 g_init = models.Gaussian1D(amplitude=a0, mean=m0, stddev=stddev0) + models.Const1D()
 fit_g = fitting.LevMarLSQFitter()
 g = fit_g(g_init, x, y)
 fwhm = 2 * g.stddev_0.value * np.sqrt(2 * np.log(2))
 return fwhm,g

def fitGauss1D_Linear(x,y,a0,m0,stddev0):
 g_init = models.Gaussian1D(amplitude=a0, mean=m0, stddev=stddev0) + models.Linear1D()
 fit_g = fitting.LevMarLSQFitter()
 g = fit_g(g_init, x, y)
 fwhm = 2 * g.stddev_0.value * np.sqrt(2 * np.log(2))
 return fwhm,g

def fitMoffat1D_Linear(x,y,a0,m0,stddev0):
 g_init = models.Moffat1D(a0,m0,stddev0) + models.Linear1D()
 fit_g = fitting.LevMarLSQFitter()
 g = fit_g(g_init, x, y)
 fwhm = 2 * g.gamma_0.value * np.sqrt(2**(1/g.alpha_0.value) - 1)
 return fwhm, g

def fitMoffat1D_Const(x,y,a0,m0,stddev0):
 g_init = models.Moffat1D(a0,m0,stddev0) + models.Const1D()
 fit_g = fitting.LevMarLSQFitter()
 g = fit_g(g_init, x, y)
 fwhm = 2 * g.gamma_0.value * np.sqrt(2**(1/g.alpha_0.value) - 1)
 return fwhm,g
 
fitmodel_dict      = {'Gauss'        : fitGauss1D,
                      'Gauss+const'  : fitGauss1D_Const,
                      'Gauss+linear' : fitGauss1D_Linear,
                      'Moffat+const' : fitMoffat1D_Const,
                      'Moffat+linear': fitMoffat1D_Linear}
fitmodel_default = 'Moffat+const'

def hfwhm(lim1, lim2, x, y):
    y_m = np.ma.masked_outside(y,lim1,lim2) 
    x_m = np.ma.masked_array(x,y_m.mask)  
    rad = np.median(np.ma.compressed(x_m))
    return rad

def profile(img, xbeg, ybeg, factor, center_method):
#  Medain filtering of noise with threshold
     fimg = ndimage.filters.median_filter(img, 3)
     dimg2 = (fimg - img)**2
     img_masked = np.ma.masked_array(img, dimg2 < fimg*factor)
     img_clean  = (img_masked - img_masked + fimg).data
     img_max    = img_clean.max()
#     printd('img_max, img.max()', img_max, img.max())
     flsort = np.sort(img_clean.ravel())
     PXSTAR = 300
     PXBKGR = 300
     pxstar = -1*PXSTAR if PXSTAR < flsort.size else -1*(flsort.size-1)
     pxbkgr = PXBKGR if PXBKGR < flsort.size else flsort.size-1
     threshold_star = flsort[pxstar]
     threshold_bkg  = flsort[pxbkgr]
# star_img = np.ma.masked_less(img_clean,threshold_star)
     nonstar_img = np.ma.masked_greater(img_clean, threshold_star)
     star_img = (nonstar_img - nonstar_img).data
     bkg_img = np.ma.masked_greater(img_clean, threshold_bkg)
     bkg = np.median(np.ma.compressed(bkg_img))

#     print("threshold_star,threshold_bkg,background =",threshold_star,threshold_bkg,bkg)
     hampl = (img_max - bkg)/2

#     print("hampl =", hampl)
     if (center_method == 'barycenter'):
      com = ndimage.measurements.center_of_mass(star_img)
#      printd('com_c_m =', com)
     elif (center_method == 'maximum_count'): 
      com = np.unravel_index(np.argmax(img), img.shape)
     else:
      com = np.unravel_index(np.argmax(img_clean),img_clean.shape)
#     printd('---------- com =', com)
#     image_centre_coord = ybeg_float+com[1], xbeg_float+com[0]
     image_centre_coord = ybeg+com[1], xbeg+com[0]    #ybeg, xbeg - python's, image_centr

     i0,j0 = com
     # img = img - threshold_bkg

     y = np.arange(img.shape[0]) - com[0]
     x = np.arange(img.shape[1]) - com[1]
     xx, yy = np.meshgrid(x, y)
     r = np.sqrt(xx**2 + yy**2).ravel()
     fl = img.ravel()
     lim1 = hampl - 0.06 * hampl
     lim2 = hampl + 0.06 * hampl
#     printd('lim1, lim2, 0.06 * hampl =', lim1, lim2, 0.06 * hampl)
     hfwhm_empir = hfwhm(lim1, lim2, r, fl-bkg)
     fwhm_empir = hfwhm_empir * 2.0
#     print('fwhm_empir =', 2.*hfwhm_empir)

     return image_centre_coord, fwhm_empir, com, r, fl
