Pardon the LaTeX, this is old documentation (and old code). You may be better off just skipping to the references section.

## The Math

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

## References

* http://en.wikipedia.org/wiki/Fractional_Fourier_transform
* D. H. Bailey and P. N. Swarztrauber, "The fractional Fourier transform and applications," SIAM Review 33, 389-404 (1991)
* http://www.gnu.org/software/gsl/manual/html_node/Fast-Fourier-Transforms.html