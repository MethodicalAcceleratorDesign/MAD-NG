.. index::
   single: elementary mathematical and physical constants

********************
Elementary Constants
********************

This chapter describes basic mathematical and physiscal constants provided by the module :code:`constant`.

Mathematical Constants
======================

This section describes basic mathematical constants uniquely defined as macros in the C header :file:`mad_cst.h` and available from C and MAD modules. If these mathematical constants are already provided by the system libraries, they are used instead of the local definitions.

==================  =====================  =========================  ======================
MAD constants       C macros               C constants                Values
==================  =====================  =========================  ======================
:code:`eps`         :macro:`DBL_EPSILON`   :const:`mad_cst_EPS`       Smallest representable increment near one
:code:`tiny`        :macro:`DBL_MIN`       :const:`mad_cst_TINY`      Smallest representable number
:code:`huge`        :macro:`DBL_MAX`       :const:`mad_cst_HUGE`      Largest representable number
:code:`inf`         :macro:`INFINITY`      :const:`mad_cst_INF`       Positive infinity, :math:`1/0`
:code:`nan`         -                      -                          Canonical NaN, :math:`0/0`
:code:`e`           :macro:`M_E`           :const:`mad_cst_E`         :math:`e, \exp(1)`
:code:`log2e`       :macro:`M_LOG2E`       :const:`mad_cst_LOG2E`     :math:`\log_2(e)`
:code:`log10e`      :macro:`M_LOG10E`      :const:`mad_cst_LOG10E`    :math:`\log_{10}(e)`
:code:`ln2`         :macro:`M_LN2`         :const:`mad_cst_LN2`       :math:`\ln(2)`
:code:`ln10`        :macro:`M_LN10`        :const:`mad_cst_LN10`      :math:`\ln(10)`
:code:`lnpi`        :macro:`M_LNPI`        :const:`mad_cst_LNPI`      :math:`\ln(\pi)`
:code:`pi`          :macro:`M_PI`          :const:`mad_cst_PI`        :math:`\pi`
:code:`twopi`       :macro:`M_2PI`         :const:`mad_cst_2PI`       :math:`2\pi`
:code:`pi_2`        :macro:`M_PI_2`        :const:`mad_cst_PI_2`      :math:`\pi/2`
:code:`pi_4`        :macro:`M_PI_4`        :const:`mad_cst_PI_4`      :math:`\pi/4`
:code:`one_pi`      :macro:`M_1_PI`        :const:`mad_cst_1_PI`      :math:`1/\pi`
:code:`two_pi`      :macro:`M_2_PI`        :const:`mad_cst_2_PI`      :math:`2/\pi`
:code:`sqrt2`       :macro:`M_SQRT2`       :const:`mad_cst_SQRT2`     :math:`\sqrt 2`
:code:`sqrt3`       :macro:`M_SQRT3`       :const:`mad_cst_SQRT3`     :math:`\sqrt 3`
:code:`sqrtpi`      :macro:`M_SQRTPI`      :const:`mad_cst_SQRTPI`    :math:`\sqrt{\pi}`
:code:`sqrt1_2`     :macro:`M_SQRT1_2`     :const:`mad_cst_SQRT1_2`   :math:`\sqrt{1/2}`
:code:`sqrt1_3`     :macro:`M_SQRT1_3`     :const:`mad_cst_SQRT1_3`   :math:`\sqrt{1/3}`
:code:`one_sqrtpi`  :macro:`M_1_SQRTPI`    :const:`mad_cst_1_SQRTPI`  :math:`1/\sqrt{\pi}`
:code:`two_sqrtpi`  :macro:`M_2_SQRTPI`    :const:`mad_cst_2_SQRTPI`  :math:`2/\sqrt{\pi}`
:code:`raddeg`      :macro:`M_RADDEG`      :const:`mad_cst_RADDEG`    :math:`180/\pi`
:code:`degrad`      :macro:`M_DEGRAD`      :const:`mad_cst_DEGRAD`    :math:`\pi/180`
==================  =====================  =========================  ======================

.. index::
   mathematical constants

Physical Constants
==================

This section describes basic physical constants uniquely defined as macros in the C header :file:`mad_cst.h` and available from C and MAD modules.

==================  =====================  =========================  ======================
MAD constants       C macros               C constants                Values
==================  =====================  =========================  ======================
:code:`minlen`      :macro:`P_MINLEN`      :const:`mad_cst_MINLEN`    Minimum length tolerance, :math:`10^-10` in :unit:`[m]`
:code:`minang`      :macro:`P_MINANG`      :const:`mad_cst_MINANG`    Minimum angle tolerance, :math:`10^-10` in :unit:`[m^{-1}]`
:code:`minstr`      :macro:`P_MINSTR`      :const:`mad_cst_MINSTR`    Minimum strength tolerance, :math:`10^-10` in :unit:`[rad]`
==================  =====================  =========================  ======================

The following table lists some physical constants from the `CODATA 2018 <https://physics.nist.gov/cuu/pdf/wall_2018.pdf>`_ sheet.

==================  =====================  =========================  ======================
MAD constants       C macros               C constants                Values
==================  =====================  =========================  ======================
:code:`clight`      :macro:`P_CLIGHT`      :const:`mad_cst_CLIGHT`    Speed of light, :math:`c` in :unit:`[m/s]`
:code:`mu0`         :macro:`P_MU0`         :const:`mad_cst_MU0`       Permeability of vacuum, :math:`\mu_0` in :unit:`[T.m/A]`
:code:`epsilon0`    :macro:`P_EPSILON0`    :const:`mad_cst_EPSILON0`  Permittivity of vacuum, :math:`\epsilon_0` in :unit:`[F/m]`
:code:`qelect`      :macro:`P_QELECT`      :const:`mad_cst_QELECT`    Elementary electric charge, :math:`e` in :unit:`[C]`
:code:`hbar`        :macro:`P_HBAR`        :const:`mad_cst_HBAR`      Reduced Plack's constant, :math:`\hbar` in :unit:`[GeV.s]`
:code:`amass`       :macro:`P_AMASS`       :const:`mad_cst_AMASS`     Unified atomic mass, :math:`m_u c^2` in :unit:`[GeV]`
:code:`emass`       :macro:`P_EMASS`       :const:`mad_cst_EMASS`     Electron mass, :math:`m_e c^2` in :unit:`[GeV]`
:code:`pmass`       :macro:`P_PMASS`       :const:`mad_cst_PMASS`     Proton mass, :math:`m_p c^2` in :unit:`[GeV]`
:code:`nmass`       :macro:`P_NMASS`       :const:`mad_cst_NMASS`     Neutron mass, :math:`m_n c^2` in :unit:`[GeV]`
:code:`mumass`      :macro:`P_MUMASS`      :const:`mad_cst_MUMASS`    Muon mass, :math:`m_{\mu} c^2` in :unit:`[GeV]`
:code:`deumass`     :macro:`P_DEUMASS`     :const:`mad_cst_DEUMASS`   Deuteron mass, :math:`m_d c^2` in :unit:`[GeV]`
:code:`eradius`     :macro:`P_ERADIUS`     :const:`mad_cst_ERADIUS`   Classical electron radius, :math:`r_e` in :unit:`[m]`
:code:`alphaem`     :macro:`P_ALPHAEM`     :const:`mad_cst_ALPHAEM`   Fine-structure constant, :math:`\alpha`
==================  =====================  =========================  ======================

.. index::
   physical constants
   CODATA
