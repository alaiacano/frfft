'''
Created on Oct 12, 2009

@author: Adam Laiacano
'''
import numpy

def nextpow2(num):
    '''
    Returns the next highest power of 2 from the given value.


    Example
    -------

    >>>nextpow2(1000)
    1024

    >>nextpow2(1024)
    2048
    '''
    npow = 2
    while npow <= num: 
        npow = npow * 2
    return npow



class FrFFT(object):
    '''
    Computes the fractional fast Fourier transform of a given input.  
    
    FrFFT(x, alpha)
    
    USAGE:
    
    import frfft
    from numpy import linspace, pi, sin
    t=linspace(-4*pi,4*pi, 1000)
    x=sin(2*pi*40*t) + sin(2*pi*20*t) + sin(2*pi*10*t)
    X=frfft.FrFFT(x,1./1024)

    subplot(211)
    plot(numpy.abs(X.result))
    subplot(212)
    plot(numpy.unwrap(arctan2(X.result.imag,X.result.real)))
    show()
    
    '''


    def __init__(self, x_in, alpha):
        '''
        Sets up required variables and computes FrFFT.
        
        @param x Input sequence
        @param alpha Fractional FFT value.
        '''
        
        ## M = length of the input vector
        self.m_in = len(x_in)
        
        ## Alpha value
        self.alpha = float(alpha * 1.)
        
        ## pp >= 2*m+1
        self.p_p = nextpow2(self.m_in + 1)
        
        ## Complex valued input data
        self.data = numpy.array(x_in, dtype = 'complex')
        
        ## Output place holder
        self.result = numpy.zeros(self.m_in, dtype = 'complex')
        
        # Run the process
        self.replace_alpha(1. * alpha)
        
        
        
    def replace_alpha(self, alpha):
        '''
        Resets the alpha value and re-runs the FrFFT.
        
        '''
        self.alpha = float(alpha)
        
        # take 2p-point FFTs of y and z
        fft_y = numpy.fft.fft(self.__generate_y(), self.p_p)
        
        fft_z = numpy.fft.fft(self.__generate_z(), self.p_p)
        
        j_g = numpy.arange(0, self.m_in, dtype = 'float')
        exponent = -1j * numpy.pi * j_g**2 * self.alpha
        y_z = numpy.fft.ifft(fft_y * fft_z, self.p_p)
        self.result = numpy.exp(exponent) * y_z[:self.m_in]
        
        
    def __generate_y(self):
        '''
        Intermediary function to calculate y described above.
        '''
        y_out = numpy.zeros(self.p_p, dtype = 'complex')
        
        j_y = numpy.arange(0, self.m_in, dtype = 'float')
        
        exponent = -1j * numpy.pi * j_y**2 * self.alpha
        
        y_out[:self.m_in] = self.data * numpy.exp(exponent)
        return y_out
    
    def __generate_z(self):
        
        '''
        Intermediary function to calculate z described above.
        '''
        z_out = numpy.zeros(self.p_p, dtype = 'complex')
        
        j_z = numpy.arange(0, self.p_p, dtype = 'float')
        pp_m = self.p_p - self.m_in        
        
        exponent = 1j * numpy.pi * j_z[:self.m_in]**2 * self.alpha
        z_out[:self.m_in] = numpy.exp(exponent)
        z_out[pp_m:] = numpy.exp(numpy.pi * j_z[pp_m:]**2 * self.alpha * 1j)

        return z_out
        
