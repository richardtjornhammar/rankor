"""
Copyright 2020 RICHARD TJÖRNHAMMAR

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
import pandas as pd
import numpy as np

def contrasted_signal ( edf ) :
    # NOT NECESSARILY VALID
    B       = np.mean ( np.mean( edf ) )
    l2_cont = lambda x,B : np.abs( ( 2**x - 2**B )/( 2**x + 2**B ) ) * x
    cdf = edf .apply( lambda x : l2_cont(x,np.mean(x)) )
    return ( cdf )

# REGULAR CONTRAST
contrast   = lambda A,B : ( A-B )/( A+B )
e_flatness = lambda x   : np.exp(np.mean(np.log(x),0))/np.mean(x,0)
e_contrast = lambda x   : 1 - e_flatness(x)

