from ...imports import *

__all__ = ["flag_outliers"]


def flag_outliers(
    self, sigma = 3
):
    """
    Flag outliers on a chromatic rainbow.

    Parameters
    ----------
    sigma : integer
        Standard deviations (sigmas) allowed for individual data points before they are flagged as outliers.

    Returns
    -------
    rainbow : Rainbow
        A new Rainbow object with the outliers flagged as 0's in the 'ok' array.
    """

    
    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("flag_outliers", locals())
    
    # create a copy of the existing rainbow
    new = self._create_copy()
    
    list_of_rows = []
    for i, row in enumerate(new.flux):
        diffs = np.diff(row)
        mask = sigma_clip(diffs, sigma=sigma, maxiters=5, cenfunc = 'median')
        
        ls = np.hstack((np.asarray([True]),~mask.mask))
        new_ls = ls.copy()
        for i in range(1,len(ls)):
            if ((ls[i] == True) and (ls[i-1] == False)):
                new_ls[i-1] = True
                
        list_of_rows.append(new_ls)
            
    full_mask = np.array(list_of_rows)

    new.ok = new.ok * full_mask  
        
    # append the history entry to the new Rainbow
    new._record_history_entry(h)
    
    return new