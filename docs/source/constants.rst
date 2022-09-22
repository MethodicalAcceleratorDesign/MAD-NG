.. index::
   Constants

*********
Constants
*********

This chapter describes some constants uniquely defined as macros in the C header :file:`mad_cst.h` and available from modules :mod:`MAD.constant` and :mod:`MAD._C` (C API) as floating point double precision variables. 

Numerical Constants
===================

.. index::
   Numerical constants

These numerical constants are provided by the system libraries. Note that the constant :const:`huge` differs from :const:`math.huge`, which corresponds in fact to :const:`inf`. 

===================  ======================  =========================  ======================
MAD constants        C macros                C constants                Values
===================  ======================  =========================  ======================
:const:`eps`         :c:macro:`DBL_EPSILON`  :const:`mad_cst_EPS`       Smallest representable step near one
:const:`tiny`        :c:macro:`DBL_MIN`      :const:`mad_cst_TINY`      Smallest representable number
:const:`huge`        :c:macro:`DBL_MAX`      :const:`mad_cst_HUGE`      Largest representable number
:const:`inf`         :c:macro:`INFINITY`     :const:`mad_cst_INF`       Positive infinity, :math:`1/0`
:const:`nan`         :c:macro:`NAN`          :const:`mad_cst_NAN`       Canonical NaN [#f1]_, :math:`0/0`
===================  ======================  =========================  ======================

Mathematical Constants
======================

.. index::
   Mathematical constants

This section describes some mathematical constants uniquely defined as macros in the C header :file:`mad_cst.h` and available from C and MAD modules as floating point double precision variables. If these mathematical constants are already provided by the system libraries, they will be used instead of their local definitions.

===================  ======================  =========================  ======================
MAD constants        C macros                C constants                Values
===================  ======================  =========================  ======================
:const:`e`           :c:macro:`M_E`          :const:`mad_cst_E`         :math:`e`
:const:`log2e`       :c:macro:`M_LOG2E`      :const:`mad_cst_LOG2E`     :math:`\log_2(e)`
:const:`log10e`      :c:macro:`M_LOG10E`     :const:`mad_cst_LOG10E`    :math:`\log_{10}(e)`
:const:`ln2`         :c:macro:`M_LN2`        :const:`mad_cst_LN2`       :math:`\ln(2)`
:const:`ln10`        :c:macro:`M_LN10`       :const:`mad_cst_LN10`      :math:`\ln(10)`
:const:`lnpi`        :c:macro:`M_LNPI`       :const:`mad_cst_LNPI`      :math:`\ln(\pi)`
:const:`pi`          :c:macro:`M_PI`         :const:`mad_cst_PI`        :math:`\pi`
:const:`twopi`       :c:macro:`M_2PI`        :const:`mad_cst_2PI`       :math:`2\pi`
:const:`pi_2`        :c:macro:`M_PI_2`       :const:`mad_cst_PI_2`      :math:`\pi/2`
:const:`pi_4`        :c:macro:`M_PI_4`       :const:`mad_cst_PI_4`      :math:`\pi/4`
:const:`one_pi`      :c:macro:`M_1_PI`       :const:`mad_cst_1_PI`      :math:`1/\pi`
:const:`two_pi`      :c:macro:`M_2_PI`       :const:`mad_cst_2_PI`      :math:`2/\pi`
:const:`sqrt2`       :c:macro:`M_SQRT2`      :const:`mad_cst_SQRT2`     :math:`\sqrt 2`
:const:`sqrt3`       :c:macro:`M_SQRT3`      :const:`mad_cst_SQRT3`     :math:`\sqrt 3`
:const:`sqrtpi`      :c:macro:`M_SQRTPI`     :const:`mad_cst_SQRTPI`    :math:`\sqrt{\pi}`
:const:`sqrt1_2`     :c:macro:`M_SQRT1_2`    :const:`mad_cst_SQRT1_2`   :math:`\sqrt{1/2}`
:const:`sqrt1_3`     :c:macro:`M_SQRT1_3`    :const:`mad_cst_SQRT1_3`   :math:`\sqrt{1/3}`
:const:`one_sqrtpi`  :c:macro:`M_1_SQRTPI`   :const:`mad_cst_1_SQRTPI`  :math:`1/\sqrt{\pi}`
:const:`two_sqrtpi`  :c:macro:`M_2_SQRTPI`   :const:`mad_cst_2_SQRTPI`  :math:`2/\sqrt{\pi}`
:const:`rad2deg`     :c:macro:`M_RAD2DEG`    :const:`mad_cst_RAD2DEG`   :math:`180/\pi`
:const:`deg2rad`     :c:macro:`M_DEG2RAD`    :const:`mad_cst_DEG2RAD`   :math:`\pi/180`
===================  ======================  =========================  ======================

Physical Constants
==================

.. index::
   Physical constants
   CODATA

This section describes some physical constants uniquely defined as macros in the C header :file:`mad_cst.h` and available from C and MAD modules as floating point double precision variables.

===============  ===================  =======================  ======================
MAD constants    C macros             C constants              Values
===============  ===================  =======================  ======================
:const:`minlen`  :c:macro:`P_MINLEN`  :const:`mad_cst_MINLEN`  Min length tolerance, default :math:`10^{-10}` in :unit:`[m]`
:const:`minang`  :c:macro:`P_MINANG`  :const:`mad_cst_MINANG`  Min angle tolerance, default :math:`10^{-10}` in :unit:`[1/m]`
:const:`minstr`  :c:macro:`P_MINSTR`  :const:`mad_cst_MINSTR`  Min strength tolerance, default :math:`10^{-10}` in :unit:`[rad]`
===============  ===================  =======================  ======================

The following table lists some physical constants from the `CODATA 2018 <https://physics.nist.gov/cuu/pdf/wall_2018.pdf>`_ sheet.

=================  =====================  =========================  ======================
MAD constants      C macros               C constants                Values
=================  =====================  =========================  ======================
:const:`clight`    :c:macro:`P_CLIGHT`    :const:`mad_cst_CLIGHT`    Speed of light, :math:`c` in :unit:`[m/s]`
:const:`mu0`       :c:macro:`P_MU0`       :const:`mad_cst_MU0`       Permeability of vacuum, :math:`\mu_0` in :unit:`[T.m/A]`
:const:`epsilon0`  :c:macro:`P_EPSILON0`  :const:`mad_cst_EPSILON0`  Permittivity of vacuum, :math:`\epsilon_0` in :unit:`[F/m]`
:const:`qelect`    :c:macro:`P_QELECT`    :const:`mad_cst_QELECT`    Elementary electric charge, :math:`e` in :unit:`[C]`
:const:`hbar`      :c:macro:`P_HBAR`      :const:`mad_cst_HBAR`      Reduced Plack's constant, :math:`\hbar` in :unit:`[GeV.s]`
:const:`amass`     :c:macro:`P_AMASS`     :const:`mad_cst_AMASS`     Unified atomic mass, :math:`m_u\,c^2` in :unit:`[GeV]`
:const:`emass`     :c:macro:`P_EMASS`     :const:`mad_cst_EMASS`     Electron mass, :math:`m_e\,c^2` in :unit:`[GeV]`
:const:`pmass`     :c:macro:`P_PMASS`     :const:`mad_cst_PMASS`     Proton mass, :math:`m_p\,c^2` in :unit:`[GeV]`
:const:`nmass`     :c:macro:`P_NMASS`     :const:`mad_cst_NMASS`     Neutron mass, :math:`m_n\,c^2` in :unit:`[GeV]`
:const:`mumass`    :c:macro:`P_MUMASS`    :const:`mad_cst_MUMASS`    Muon mass, :math:`m_{\mu}\,c^2` in :unit:`[GeV]`
:const:`deumass`   :c:macro:`P_DEUMASS`   :const:`mad_cst_DEUMASS`   Deuteron mass, :math:`m_d\,c^2` in :unit:`[GeV]`
:const:`eradius`   :c:macro:`P_ERADIUS`   :const:`mad_cst_ERADIUS`   Classical electron radius, :math:`r_e` in :unit:`[m]`
:const:`alphaem`   :c:macro:`P_ALPHAEM`   :const:`mad_cst_ALPHAEM`   Fine-structure constant, :math:`\alpha`
=================  =====================  =========================  ======================

.. ------------------------------------------------------------

.. rubric:: Footnotes

.. [#f1] Canonical NaN bit patterns may differ between MAD and C for the mantissa, but both should exibit the same behavior.
