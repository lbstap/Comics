###########################################################################################################################
###                                COnceptual Model of Ice-sheet ChangeablenesS (COMICS)                                ###
###########################################################################################################################
# Execute: $ python comics.py
#
# Program to calculate transient ice volume of an ice sheet, depending on a prescribed
# control parameter - equilibrium ice volume (C-Veq) relation and evolution of the control parameter
# Basic procedure:
# If the ice volume is smaller than the equilibrium volume, it grows by the growth rate.
# If the ice volume is larger than the equilibrium volume, it shrinks by the decay rate.
#
# Cite:
#
#  Stap, L.B., Knorr, G., and Lohmann, G.: Anti-phased Miocene ice volume and CO2 changes by transient Antarctic ice sheet
#    variability, Paleoceanography and Paleoclimatology, in review.
#
# -------------------------------------------------------------------------------------------------------------------------
# Written by: L.B. Stap at Alfred-Wegener-Institut, Helmholtz-Zentrum fuer Polar- und Meeresforschung, Bremerhaven, Germany
# Email: lennert.stap<at>gmail.com
# July 2019, last change: October 2020
