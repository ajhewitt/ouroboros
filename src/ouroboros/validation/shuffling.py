"""
Shuffling Engine (Plan C Validation)
Generates randomized catalogs to rule out Window Function aliasing.
"""

import numpy as np

def shuffle_catalog_vectors(ra, dec, method='spin'):
    """
    Randomizes quasar positions or vectors to create a null set.
    
    Args:
        ra, dec (array): Original coordinates.
        method (str): 
            'spin': Rotates the longitude (RA) randomly. Preserves latitude distribution.
            'scramble': Randomly re-assigns positions from the dataset to different objects.
            
    Returns:
        tuple: (ra_null, dec_null)
    """
    n_objs = len(ra)
    rng = np.random.default_rng()
    
    if method == 'spin':
        # Apply a random RA offset to everything (keeping Dec fixed)
        # This preserves the "clumping" in Dec but spins the survey around the pole.
        # Good for checking longitudinal alignments like the "Axis of Evil".
        ra_offset = rng.uniform(0, 360)
        ra_null = (ra + ra_offset) % 360
        return ra_null, dec
        
    elif method == 'scramble':
        # complete reshuffle
        indices = np.arange(n_objs)
        rng.shuffle(indices)
        return ra[indices], dec[indices]
        
    return ra, dec
