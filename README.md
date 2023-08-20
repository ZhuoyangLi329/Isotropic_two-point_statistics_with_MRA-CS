# Isotropic two-point statistics with MRA-CS



#### The distribution of galaxies

let‘s start from the mathematical form of the distribution of galaxies.

$$
n(\mathbf{x})=\sum_{i=0}^N w_i \delta_D^3\left(\mathbf{x}-\mathbf{x}_{\mathbf{i}}\right)
$$

where $\delta_D$ is Dirac delta function, $N$ is the total number of sampling galaxies, $\mathbf{x}_i$ is the position of the $i$ th galaxy and $W_i$ is its weight.

#### Filtered Density Field

For most of measurable statistical quantities, we need to evaluate

$$
n_W(\mathbf{x})=\int W\left(\mathbf{x}-\mathbf{x}^{\prime}\right) n\left(\mathbf{x}^{\prime}\right) d^3 \mathbf{x}^{\prime}
$$

where $W(\mathbf{x})$ is a kernel. In the $k$-space we have 

$$
\hat{n}_W(\mathbf{k})=\hat{W}(\mathbf{k}) \hat{n}(\mathbf{k})
$$

We list some mathematical expressions for kernel in real and $k$-space.

##### Spherical Top Hat

$$
W_{\text {sphere }}(r, R)=\frac{1}{(4 \pi / 3) R^3} \theta(R-r)
$$

$$
\hat{W}_{\text {sphere }}(k, R)=\frac{3}{k^3 R^3}(\sin (k R)-k R \cos (k R))
$$

##### Spherical Shell Top Hat

$$
W_{\text {shell }}(r, R)=\frac{1}{4 \pi r^2} \delta_D(R-r)
$$

$$
W_{\text {shell }}(k, r)=\frac{\sin (k r)}{k r}
$$



##### Cubic Top Hat

$$
W_{\text {cubic }}(\mathbf{x}, \mathbf{L})=\frac{1}{L_x L_y L_z} \theta\left(L_x / 2-|x|\right) \theta\left(L_y / 2-|y|\right) \theta\left(L_z / 2-|z|\right)
$$

$$
\hat{W}_{\text {cubic }}(\mathbf{k}, \mathbf{L})=\frac{\sin \left(k_x L_x / 2\right)}{k_x L_x / 2} \frac{\sin \left(k_y L_y / 2\right)}{k_y L_y / 2} \frac{\sin \left(k_z L_z / 2\right)}{k_z L_z / 2}
$$

where the vector $\mathbf{L}$ denoting for ${L_x, L_y, L_z}$ defines the side length of cubic.

##### Cylinder Top Hat

$$
W_{\text {cylinder }}(\rho, z, R, h)=\frac{1}{\pi R^2 h} \theta(R-\rho) \theta(h / 2-|z|)
$$

$$
\hat{W}_{\text {cylinder }}\left(k_{\perp}, k_z, R, h\right)=\frac{\sin \left(k_z h / 2\right)}{k_z h / 2} \int_0^1 J_0\left(k_{\perp} R \sqrt{x}\right) d x
$$

where $k_{\perp}=\sqrt{k_x^2+k_y^2}$

##### Gaussian Filter

$$
W_G(r, R)=\frac{1}{\sqrt{(2 \pi)^3} R^3} \exp \left(-\frac{r^2}{2 R^2}\right)
$$

$$
\hat{W}_G(k, R)=\exp \left(-\frac{1}{2} k^2 R^2\right)
$$

### Multi-Resolution Analysis

#### MRA-Wavelet

When we select a resolution j, we get a series of wavelet basis functions that can be used to represent an arbitrary function.

$$
\phi_{j, k}(x)=2^{j / 2} \phi\left(2^j x-k\right) \quad \mid k \in \mathbf{Z}
$$

where $\phi$ is the basic scaling function.

#### MRA-Density Field

We can make a decomposition of the density distribution $n(x)=\sum w_i \delta\left(x-x_i\right)$ in the MRA at a scale $j$ .

$$
n(x)=\sum_l s_l^j \phi_{j, l}(x)
$$

where $s_1^j$ is the scaling function coefficients (SFCs).

$$
\begin{aligned}
s_l^j & =\int n(x) \phi_{j, l}(x) d x \\
& =\sum_{i=1}^{N_p} w_i 2^{j / 2} \phi\left(2^j x_i-I\right)
\end{aligned}
$$

#### MRA-Kernel

We can make a decomposition of the kernel $W(x, y)$ (a 2D example is given here)

$$
W_j(x, y)=\sum_{l, m} w_{l, m}^j \phi_{j, I}(x) \phi_{j, m}(y)
$$

where

$$
w_{l, m}^j=\int W(x, y) \phi_{j, l}(x) \phi_{j, m}(y) d x d y
$$

#### MRA-Filtered density field

We can make a decomposition of the filtered density field $n_W(x)$ 

$$
n_W(x) \rightarrow n_W^j(x)=\sum_I \tilde{s}_l^j \phi_{j, l}(x)
$$

with

$$
\tilde{\boldsymbol{s}}_l^j=\sum_m w_{l, m}^j \boldsymbol{s}_m^{\dot{j}}
$$

For a homogenous kernel $W , w_{l, m}^j=w_{l-m}^j$ is a Toeplitz matrix, and conventionally its operation on the vector $\mathbf{s}^j=\left\{\boldsymbol{s}_m^j\right\}$ could be accomplished using the fast Fourier transformation technique (FFT) **Here we used the Beylkin & Cramer Summation Rule**

### MRACS-One point statistics

 One-point statistics are not the focus of this work, so we only list here the results that have been available previously 

##### One point statistics- Count-In-Cell (CIC)

![Fig_cic](C:\Users\lzy\Downloads\Fig_cic.png)

##### One point statistics- PDF

![1692363526338](C:\Users\lzy\Documents\WeChat Files\wxid_g9vkpgy4v9ws21\FileStorage\Temp\1692363526338.png)

### MRACS-Two point statistics

#### Two point statistics-$\sigma$

The second order variance of filtered density fluctuations $\sigma$ can be calculated by

$$
\sigma^2=\frac{1}{(2 \pi)^3} \int\left|W_{\text {window }}(\mathbf{k})\right|^2 P(k) d^3 \mathbf{k}
$$

If we choose spherical kernel $\hat{W}_{\text {sphere }}(k, R)=\frac{3}{k^3 R^3}(\sin (k R)-k R \cos (k R))$, we will return to the traditional definition of $\sigma$ .

Using the definition of filtered density field, we can use wavelet coefficients to represent $\sigma$ 

$$
\sigma^2(\cdot)=\sum_l \tilde{\mathbf{s}}_l^j \tilde{\mathbf{s}}_l^j=\left|\tilde{\mathbf{s}}^j\right|^2-1
$$

Below, we present the theoretical results of traditional $\sigma$ and the calculation results of MRA-CS in redshift space. (We use the 

[nbodykit]: https://nbodykit.readthedocs.io/en/latest/index.html

 software package to calculate the theoretical Halofit nonlinear power spectrum)**(simulation:*MDPL2, cdm*)**

![Fig_tra_sigma](C:\Users\lzy\Downloads\Fig_tra_sigma.jpg)

Next, we choose the cylindrical window to calculate the variance.

$$
\hat{W}_{\text {cylinder }}\left(k_{\perp}, k_z, R, h\right)=\frac{\sin \left(k_z h / 2\right)}{k_z h / 2} \frac{2 J_1\left(k_{\perp} R\right)}{k_{\perp} R}
$$

where $k_{\perp}=\sqrt{k_x^2+k_y^2}$

In real space we have

$$
\begin{aligned}
& \sigma_{\text {real-space }}^2=\frac{1}{(2 \pi)^3} \int\left|\frac{\sin \left(k_z h / 2\right)}{k_z h / 2} \frac{2 J_1(k R)}{k R}\right|^2 P_{\text {nonlin }}(k) d^3 \mathbf{k} \\
& =\frac{1}{(2 \pi)^2} \int\left(\frac{\sin \left(k_z h / 2\right)}{k_z h / 2}\right)^2 d k_z \int P_{\text {nonlin }}(k)\left(\frac{2 J_1\left(k_{\perp} R\right)}{k_{\perp} R}\right)^2 k_{\perp} d k_{\perp} \\
&
\end{aligned}
$$

In redshift space we have

$$
\begin{aligned}
& \sigma_{\text {redshift-space }}^2=\frac{1}{(2 \pi)^3} \int\left|\frac{\sin \left(k_z h / 2\right)}{k_z h / 2} \frac{2 J_1(k R)}{k R}\right|^2 P_{\text {nonlin }}^s(k, \mu) d^3 \mathbf{k} \\
= & \frac{1}{(2 \pi)^2} \int\left(\frac{\sin \left(k_z h / 2\right)}{k_z h / 2}\right)^2 d k_z \int P_{\text {nonlin }}^s(k, \mu)\left(\frac{2 J_1\left(k_{\perp} R\right)}{k_{\perp} R}\right)^2 k_{\perp} d k_{\perp}
\end{aligned}
$$

where $P^{\mathrm{s}}(k, \mu)=\left(b+f \mu^2\right)^2 P_{\text {nonlin }}(k) D^{\mathrm{FoG}}\left(k \mu \sigma_v\right)$ represents redshift-space galaxy power spectrum, $\mu=k_z / k$ and $D^{\mathrm{FoG}}=\left(1+k^2 \mu^2 \sigma_v^2 / H^2\right)^{-1}$

The following figure shows an example in a redshift space: comparison between theoretical calculations and MRA-CS calculations (in contour form). The R and H in the figure represent the radius and height of the cylindrical kernel, respectively **(simulation:*Quijote, halo*)**

![9867a66b0ce6789c3127136bfd88aeb](C:\Users\lzy\AppData\Local\Temp\WeChat Files\9867a66b0ce6789c3127136bfd88aeb.jpg)

##### $\sigma$-parameter estimation

For example, we choose 10 data points $\left(r_i, h_i\right), i=0-9$ and calculate $\sigma$ or $d \sigma / d h$ for these data points with MRA-CS in 1000 halo samples, now we have 1000 data vectors $\boldsymbol{\xi}^{\text {data }}$

we calculate the covariance matrix of 1000 halo samples.

$$
\mathbf{C}=\frac{1}{N-1} \sum_{k=1}^N\left(\boldsymbol{\xi}_k^{\text {data }}-\overline{\boldsymbol{\xi}^{\text {data }}}\right)\left(\boldsymbol{\xi}_k^{\text {data }}-\overline{\boldsymbol{\xi}^{\text {data }}}\right)
$$

where $N=1000$, and $\overline{\boldsymbol{\xi}^{\text {data }}}$ is the mean data vector.

The inversion of covariance matrix provides a biased estimator of the true precision matrix. In order to account for this effect, we multiply $\mathbf{C}^{-1}$ by a correction factor $\alpha=\left(1-\frac{N_{\mathrm{b}}+1}{N-1}\right)$ ***(Hartlap et al. 2007)***, where $N_b=10, N=1000$.

The $\chi^2$ can be calculated by

$$
\chi^2=\left(\boldsymbol{\xi}^{\text {theory }}-\boldsymbol{\xi}^{\text {data }}\right) \mathbf{C}^{-1}\left(\boldsymbol{\xi}^{\text {theory }}-\boldsymbol{\xi}^{\text {data }}\right)^T
$$

Then, the log likelihood can be expressed as

$$
\log \mathcal{L}=-\frac{1}{2} \chi^2
$$

After obtaining the likelihood function, we use the Python software package 

[emcee]: https://emcee.readthedocs.io/en/stable/

 for parameter estimation.

However, we encountered a failure here, and the following figure shows us why $\sigma$ is not suitable for parameter estimation (The correlation between different data points is too high).

![Fig_sigma_faild](C:\Users\lzy\Downloads\Fig_sigma_faild.jpg)

#### Two point statistics-correlation function $\xi$

#### 2PCF

Let's start with the simplest correlation function -2 point correlation function.

The two-point correlation function $\xi(\mathbf{r})=<\delta(\mathbf{x}) \delta(\mathbf{x}+\mathbf{r})>$ is Fourier conjugate of the power spectrum $\left.P(k)=<|\delta(\mathbf{k})|^2\right\rangle$.

$$
\xi(r)=\frac{1}{2 \pi^2} \int_0^{\infty} \hat{W}_{\text {shell }}(k, r) P(k) k^2 d k
$$

We use the LS estimator to calculate the 2-point correlation function.

$$
\hat{\xi}_{L S}(r)=\frac{D D-2 D R+R R}{R R}
$$

$D D$ : the count of pairs in the catalog within the interval $[r-0.5 d r, r+0.5 d r]$

$R R$ : the number of pairs in a random sample in the same interval; 

$D R$ : the number of pairs between the catalog and the binomial random sample in the same interval.

The pair count in the shell of $(r, r+d r)$ from a particle $i$ is approximated by

$$
n_W^j\left(x_i\right)=d V \sum_I \tilde{s}_l^j \phi_{j, l}\left(x_i\right)
$$

where $d V=4 \pi r^2 d r$ is the differential volume element of spherical shell if $d r$ is small enough
Summing over the whole sample, the total number of pairs is

$$
D D=\sum_{i=1}^N w_i \sum_l \tilde{s}_l^j \phi_{j, l}\left(x_i\right) \Longrightarrow D D=\sum_l \tilde{s}_l^j \int d x \phi_{j, l}(x) n(x)
$$

where $n(x)$ is the number density. Using the orthogonality of the base functions, we have

$$
D D=\sum_l \tilde{s}_l^j s_l^j=\tilde{\mathbf{s}}^j \cdot \mathbf{s}^j
$$

$DR$ and $RR$ can be calculated using the same method.

The following figure shows the comparison of 2PCF calculated by MRA-CS with the official 2PCF provided by "Quijote Simulation"（Note that the horizontal axis of the graph is $r(Mpc/h)$, the vertical axis is $r^2\xi$）

![Fig_qjt_2pcf_real](C:\Users\lzy\Downloads\Fig_qjt_2pcf_real.jpg)

#### 2D-2PCF

In redshift space, the 2PCF has an additional angular parameter, so the 2D-2PCF can be viewed as an correlation of particles at a specific distance and angle, i.e., a "ring" on the spherical shell.

$$
\hat{W}_{\text {ring }}(k, R,\alpha)=J_0\left(k_{\perp} R \sin \alpha\right) e^{-i k_z R \cos \alpha}
$$

where $\alpha$ is the angle between the particle direction and the line of sight direction (z-direction).

The RSD effect destroys the symmetry of the particle distribution, which is manifested in the 2D-2PCF, where the ring kernel loses one direction of symmetry compared to the spherical shell kernel. This also predetermines the superiority of 2D-2PCF over conventional 2PCF in extracting cosmological information in redshift space.

To compute 2D-2PCF with MRA-CS, simply replace the spherical shell kernel with the ring kernel, then perform the same steps as for 2PCF.

In calculating the 2D-2PCF through theory, we refer to the 

[1511.00012]: https://arxiv.org/abs/1511.00012

(eq.26-eq.32)

#### "Cylinder/Disk" correlation function

Since the distribution of particles is non-isotropic in redshift space, choosing a non-isotropic kernel for the correlation of particles will be more helpful for us to extract cosmological information. MRA-CS can help us to do this very easily (it supports arbitrary shaped kernels, and it doesn't even need a resolved form).

In this work, we are not going to use a mathematically very complex kernel, but rather to illustrate the superiority of MRA-CS through some easy-to-understand examples, outside of the 2D-2PCF (ring kernel), where we consider that particles can also be associated with other particles through the side or top faces of a cylinder.

We first find the form of the cylindrical side surface (called "cylshell" in the following) and the cylindrical top surface (called "disk" in the following) in the real space and $k$-space.

$$
W_{\text {cylshell }}(r, z, R, H)=\frac{1}{4 \pi R H} \delta_D(R-r) \theta(H-|z|)
$$

$$
\hat{W}_{\text {cylshell }}(R, H, \vec{k})=\frac{\sin \left(k_z H\right)}{k_z H} J_0(k_\perp R)
$$

$$
W_{\text {disk }}(r, z, R, H)=\frac{1}{2 \pi R^2} \theta(R-r) \delta_D(H-|z|)
$$

$$
\hat{W}_{\text {disk }}(R, H, \vec{k})=\frac{2 \cos \left(k_z H\right) J_1\left(k_{\perp} R\right)}{k_{\perp} R}
$$

Where R and H are the radius and half-height of the cylinder, respectively.

To compute these two correlation functions with MRA-CS, one simply replaces the corresponding kernel, and the rest of the steps are exactly the same as in the previous computation of the 2PCF.

The next problem we encounter is how to theoretically compute these two correlation functions. In fact, they can be calculated by integrating the 2D-2PCF over both surfaces.

$$
\xi_{\text {cylshell }}(r,\mu)=\iiint \xi{((z^2+\rho^2)^{1/2}, z/(z^2+\rho^2)^{1/2})} W_{\text {cylshell }}(\rho, z, R, H) \rho d \rho d \theta d z
$$

$$
\xi_{\text {disk }}(r,\mu)=\iiint \xi{((z^2+\rho^2)^{1/2}, z/(z^2+\rho^2)^{1/2})} W_{\text {disk }}(\rho, z, R, H) \rho d \rho d \theta d z
$$

where $\xi$ is 2D-2PCF,  $r^2=(R^2+H^2)$ is half-length of the diagonal of the cylinder, $\mu=H/\rho$ is the cosine of the angle with the direction of the line of sight.

In order to verify the consistency between the theoretical calculations and the MRA-CS calculations, we calculated the "Cylinder/Disk" correlation function in 1000 halo samples (***Quijote simulation***) with MRA-CS and averaged the results. It should also be noted that a linear RSD correction is used in the theoretical calculation of the 2D-2PCF for speed reasons.

![Fig_xi_cyl](C:\Users\lzy\Downloads\Fig_xi_cyl.jpg)

![Fig_xi_disk](C:\Users\lzy\Downloads\Fig_xi_disk.jpg)

#### Parameter estimation

Now we turn to the parameter estimation of these correlation functions, we compute the RSD-2PCF, Cylinder correlation function, Disk correlation function,and extracted their multipole moments $\ell$=0,2. Each multipole, with 20 data points, corresponds to a uniform sampling of r from 50 Mpc/h to 160 Mpc/h. For each correlation function, we used 1000 halo samples to calculate the covariance matrix. The steps for parameter estimation are consistent with those described in the $\sigma$-parameter estimation section.

Since a fast code has not yet been found for the computation of the nonlinear 2D-2PCF, we can only use the linear 2D-2PCF for the theoretical computations, which means that we will lose the nonlinear details and are currently limited to consider only the two parameters - the linear growth factor $f$ and the bias $b$.

The following three graphs show the results of parameter constraints of 2D-2PCF,Cylinder correlation function, and Disk correlation function in order.

![Fig_para_2d2pcf](C:\Users\lzy\Downloads\Fig_para_2d2pcf.jpg)

![fig_para_cyls](C:\Users\lzy\Downloads\fig_para_cyls.jpg)

![Fig_para_disk](C:\Users\lzy\Downloads\Fig_para_disk.jpg)

Next we will introduce the 2PCF nonlinear correction (with the new parameter $\sigma_v$) to confirm whether these correlation functions really contain different cosmological information.
