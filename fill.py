import numpy as np
import xarray as xr
from numba import jit


def lateral_fill(da_in, isvalid_mask, ltripole=False, tol=1.0e-4,
                 use_sor=False, rc=1.8, max_iter=1000):
    """Perform lateral fill on xarray.DataArray

    Parameters
    ----------
    da_in : xarray.DataArray
      DataArray on which to fill NaNs. Fill is performed on the two
      rightmost dimenions. Grid is assumed periodic in `x` direction
      (last dimension).

    isvalid_mask : xarray.DataArray, boolean
      Valid values mask: `True` where data should be filled. Must have the
      same rightmost dimenions as `da_in`.

    ltripole : boolean, optional [default=False]
      Logical flag; if `True` then treat the top row of the grid as periodic
      in the sense of a tripole grid.

     tol : float, optional [default=1.0e-4]
      Convergence criteria: stop filling when values change is less or equal
      to `tol * var`; i.e. `delta <= tol * np.abs(var[j, i])`.

    use_sor: boolean, optional [default=False]
      switch to select SOR fill algorithm over progressive fill algorithm

    rc : float, optional [default=1.8, valid bounds=(1.0,2.0)]
       over-relaxation coefficient to use in SOR fill algorithm. Larger arrrays
       typically converge faster with larger coefficients. For 1 deg. grid (360x180)
       a coefficient in the range 1.85-1.9 is near optimal.

    max_iter : integer, optional, [default=1000]
       maximum number of iterations to do before giving up if tol is not reached.

    Returns
    -------
    da_out : xarray.DataArray
      DataArray with NaNs filled by iterative smoothing.

    """
    print("IN FOB version : lateral_fill")

    dims_in = da_in.dims
    non_lateral_dims = dims_in[:-2]

    attrs = da_in.attrs
    encoding = da_in.encoding
    coords = da_in.coords

    da_in, isvalid_mask = xr.broadcast(da_in, isvalid_mask)

    if len(non_lateral_dims) > 0:
        da_in_stack = da_in.stack(non_lateral_dims=non_lateral_dims)
        da_out_stack = xr.full_like(da_in_stack, fill_value=np.nan)
        isvalid_mask_stack = isvalid_mask.stack(non_lateral_dims=non_lateral_dims)
        for i in range(da_in_stack.shape[-1]):
            arr = da_in_stack.data[:, :, i]
            da_out_stack[:, :, i] = lateral_fill_np_array(arr, isvalid_mask_stack.data[:, :, i],
                                                          ltripole,tol,use_sor,rc,max_iter)

        da_out = da_out_stack.unstack('non_lateral_dims').transpose(*dims_in)

    else:
        da_out = xr.full_like(da_in, fill_value=np.nan)
        da_out[:, :] = lateral_fill_np_array(da_in.data, isvalid_mask.data,
                                             ltripole,tol,use_sor,rc,max_iter)

    da_out.attrs = attrs
    da_out.encoding = encoding
    for k, da in coords.items():
        da_out[k].attrs = da.attrs

    return da_out


def lateral_fill_np_array(var, isvalid_mask, ltripole=False, tol=1.0e-4,
                          use_sor=False, rc=1.8, max_iter=1000):

    """Perform lateral fill on numpy.array

    Parameters [NB defaults set redundantly with lateral_fill above to allow this
                function to be called directly for numpy arrays]
    ----------

    var : numpy.array
      Array on which to fill NaNs. Fill is performed on the two
      rightmost dimenions. Grid is assumed periodic in `x` direction
      (last dimension). Only NaNs where isvalid_mask is True will be filled.

    isvalid_mask : numpy.array, boolean
      Valid values mask: `True` where data should be filled. Must have the
      same rightmost dimenions as `da_in`.

    ltripole : boolean, optional [default=False set in lateral_fill]
      Logical flag; if `True` then treat the top row of the grid as periodic
      in the sense of a tripole grid.

    tol : float, optional [default=1.0e-4 set in lateral_fill]
      Convergence criteria: stop filling when values change is less or equal
      to `tol * var`; i.e. `delta <= tol * np.abs(var[j, i])`.

    use_sor: boolean, optional [default=False set in lateral_fill]
      switch to select SOR fill algorithm

    rc : float, optional [default=1.8, valid bounds=(1.0,2.0) set in lateral_fill]
       over-relaxation coefficient to use in SOR fill algorithm. Larger arrrays
       typically converge faster with larger coefficients.

    max_iter : integer, optional, [default=1000 set in lateral_fill]
       maximum number of iterations to do before giving up if tol is not reached.

    Returns
    -------

    da_out : xarray.DataArray
      DataArray with NaNs filled by iterative smoothing.

    """

#    print("IN FOB version : lateral_fill_np_array")

    fillmask = np.isnan(var) & isvalid_mask
    nlat, nlon = var.shape[-2:]
    missing_value = 1e36

    var = var.copy()

    if use_sor :
        _iterative_fill_sor(nlat, nlon, var, fillmask, tol, rc, max_iter, ltripole)
    else :
        var[np.isnan(var)] = missing_value
        _iterative_fill_POP_core(nlat, nlon, var, fillmask, missing_value, tol, ltripole)

    var[var == missing_value] = np.nan

    return var


@jit(nopython=True)
def _iterative_fill_POP_core(nlat, nlon, var, fillmask, missing_value, tol, ltripole):
    """Iterative smoothing algorithm."""

#    print("IN FOB version : _iterative_fill_POP_core")
    done = False
    iter = 0

    work = np.empty((nlat, nlon))

    while not done:
        done = True
        iter += 1

        # assume bottom row is land, so skip it
        for j in range(1, nlat):
            jm1 = j - 1
            jp1 = j + 1

            for i in range(0, nlon):
                # assume periodic in x
                im1 = i - 1
                if i == 0:
                    im1 = nlon - 1
                ip1 = i + 1
                if i == nlon - 1:
                    ip1 = 0

                work[j, i] = var[j, i]

                if not fillmask[j, i]:
                    continue

                numer = 0.0
                denom = 0.0

                # East
                if var[j, ip1] != missing_value:
                    numer += var[j, ip1]
                    denom += 1.0

                # North
                if j < nlat - 1:
                    if var[jp1, i] != missing_value:
                        numer += var[jp1, i]
                        denom += 1.0

                else:
                    # assume only tripole has non-land top row
                    if ltripole:
                        if var[j, nlon - 1 - i] != missing_value:
                            numer += var[j, nlon - 1 - i]
                            denom += 1.0

                # West
                if var[j, im1] != missing_value:
                    numer += var[j, im1]
                    denom += 1.0

                # South
                if var[jm1, i] != missing_value:
                    numer += var[jm1, i]
                    denom += 1.0

                # self
                if var[j, i] != missing_value:
                    numer += denom * var[j, i]
                    denom *= 2.0

                if denom > 0.0:
                    work[j, i] = numer / denom
                    if var[j, i] == missing_value:
                        done = False
                    else:
                        delta = np.fabs(var[j, i] - work[j, i])
                        if delta > tol * np.abs(var[j, i]):
                            done = False

        var[1:nlat, :] = work[1:nlat, :]

@jit(nopython=True)
def _iterative_fill_sor(nlat, nlon, var, fillmask, tol=5.0e-4,
            rc=1.6, max_iter=100, ltripole=False):
    """Iterative land fill algorithm via SOR solution of Laplace Equation."""

#    print("IN FOB version : _iterative_fill_sor")

    # Compute a zonal mean to use as a first guess
    # Apprarently jit doesn't like masked arrays so loop it out
    zoncnt = np.zeros(nlat)
    zonavg = np.zeros(nlat)
    for j in range(0,nlat) :
        zoncnt[j] = np.sum(np.where(fillmask[j,:],0,1))
        zonavg[j] = np.sum(np.where(fillmask[j,:],0,var[j,:]))
        if zoncnt[j] != 0 : zonavg[j] = zonavg[j]/zoncnt[j]

    # Fill missing zonal averages for rows that are entirely land
    for j in range(0,nlat-1) :   # northward pass
        if zoncnt[j] > 0 and zoncnt[j+1] == 0:
            zoncnt[j+1]=1
            zonavg[j+1] = zonavg[j]
    for j in range(0,nlat-1) :  # southward pass
        jrev = nlat-1-j
        if zoncnt[jrev] > 0 and zoncnt[jrev-1] == 0 : 
            zoncnt[jrev-1]=1
            zonavg[jrev-1] = zonavg[jrev]

    # Replace the input array missing values with zonal average as first guess
    for j in range(0,nlat) :
        for i in range(0,nlon) :
            if fillmask[j,i] : var[j,i] = zonavg[j]

    # Now do the iterative 2D fill
    res = np.zeros((nlat,nlon))  # work array hold residuals
    res_max = tol
    iter = 0
    while iter < max_iter and res_max >= tol:
        res = res*0.0  # reset the residual to zero for this iteration

        # assume bottom row is all land, leave it set to zonal average
        # deal with top row separately below
        for j in range(1, nlat-1):
            jm1 = j - 1
            jp1 = j + 1

            for i in range(0, nlon):
                if fillmask[j, i]:
                    im1 = i - 1
                    if i == 0:                  # assume periodic in x
                        im1 = nlon - 1
                    ip1 = i + 1
                    if i == nlon - 1:
                        ip1 = 0

                    # this is SOR
                    res[j,i] = var[j,ip1] + var[j,im1] + var[jm1,i] + var[jp1,i] - 4.0*var[j,i]
                    var[j,i] = var[j,i] + rc*0.25*res[j,i]

        # do the top row if there was some valid data there in the input
        # otherwise leave it set to zonal average of northernmost row with valid data
        if  zoncnt[nlat-1] > 1 :
            j = nlat-1
            jm1 = j-1
            jp1 = j
            for i in range(0,nlon) :
                if fillmask[j,i] :
                    im1 = i-1
                    if i == 0:
                        im1 = nlon - 1
                    ip1 = i+1
                    if i == nlon - 1:
                        ip1 = 0
                    io = nlon-1-i

                    if ltripole :  # use cross-pole periodicity
                        res[j,i] = var[j,ip1] + var[j,im1] + var[jp1,io] + var[jm1,i] - 4.0*var[j,i]
                        var[j,i] = var[j,i] + rc*0.25*res[j,i]
                    else:          # do a 1D smooth on pole row
                        res[j,i] = var[j,ip1] + var[j,im1] - 2.0*var[j,i]
                        var[j,i] = var[j,i] + rc*0.5*res[j,i]


        res_max = np.max(np.abs(res))
        iter += 1
        
    return (iter,res_max)
