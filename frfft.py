'''
Created on Oct 12, 2009

@author: Adam Laiacano
'''
import numpy

class FrFFT(object):
    '''
    Computes the fractional fast Fourier transform of a given input.  
    
    FrFFT(x, alpha)
    
    USAGE:
    
    from gsbra.dsp import frfft
    from numpy import linspace, pi, sin
    t=linspace(-4*pi,4*pi, 1000)
    x=sin(2*pi*40*t) + sin(2*pi*20*t) + sin(2*pi*10*t)
    X=frfft.FrFFT(x,1./1024)

    subplot(211)
    plot(numpy.abs(X.result))
    subplot(212)
    plot(numpy.unwrap(arctan2(X.result.imag,X.result.real)))
    show()
    
    BACKGROUND:
    
    The FrFFT equation is \f$ G_k(\textbf{\textbf{x}},\alpha) =
    \displaystyle\sum_{j=0}^{m-1} x_{j}e^{-2 \pi i j k \alpha} \f$.
    
    To calculate it, re-write the equation as
    
    \f$ G_k(\textbf{\textbf{x}},
    \alpha) = e^{\pi i k^2 \alpha} \displaystyle\sum_{j=0}^{m-1} y_{j} z_{k-j}
    \f$
    
    where the m-long sequences y and z are defined by
    
    \f$\\y_j = x_j e^{- \pi i j^2 \alpha} \\
    z_j = e^{\pi i j^2 \alpha}
    \f$
    
    To ensure a radix-2 FFT, we extend the sequences to a length of \f$ 2p =
    nextpow2(m+1) \f$ according to the following:
    
    \f$
    \begin{tabular}{ l r }
     $y_j = 0$ & $m \le j < 2p $ \\
     $z_j = 0$ & $m \le j < 2p-m$ \\
     $z_j = e^{\pi i (j-2p)^2 \alpha}$ & $ 2p-m \le j < 2p $ \\
       \end{tabular}
    \f$
    
    This satisfies the properties for a 2p-point circular convolution and
    standard FFTs can be used.  The formula becomes
    
    \f$ G_k(\textbf{\textit{x}}, \alpha) = e^{- \pi i k^2 \alpha} F_k^{-1}
    (\textbf{\textit{w}})\f$
    
    Where \f$\textbf{w}\f$ is the 2p-long sequence defined by
    \f$\textit{w} = F_k(\textbf{y}) F_k(\textbf{z}) \f$
    
    Therefore, the Fractional FFT becomes
    
    \f$ G_k(\textbf{\textit{x}}, \alpha) = e^{- \pi i k^2 \alpha} F_k^{-1}
    \left(F_k(\textbf{y}) F_k(\textbf{z})\right) \f$
    
    I use the GNU Statistical Library's FFT functions and compared the output
    to what I calculated in Matlab.
    
    \sa http://en.wikipedia.org/wiki/Fractional_Fourier_transform
    \sa D. H. Bailey and P. N. Swarztrauber, "The fractional Fourier
    transform and applications," SIAM Review 33, 389-404 (1991)
    \sa http://www.gnu.org/software/gsl/manual/html_node/Fast-Fourier-Transforms.html
    
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
        self.p_p = self.nextpow2(self.m_in + 1)
        
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
        
        
    def generate_y(self):
        '''
        Intermediary function to calculate y described above.
        '''
        y_out = numpy.zeros(self.p_p, dtype = 'complex')
        
        j_y = numpy.arange(0, self.m_in, dtype = 'float')
        
        exponent = -1j * numpy.pi * j_y**2 * self.alpha
        
        y_out[:self.m_in] = self.data * numpy.exp(exponent)
        return y_out
    
    def generate_z(self):
        
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

