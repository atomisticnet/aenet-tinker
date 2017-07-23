c
c
c     ############################################################
c     ##      COPYRIGHT (C)  1990  by  Jay William Ponder       ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine extra1  --  user defined extra potentials  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "extra1" calculates any additional user defined potential
c     energy contribution and its first derivatives
c
c ----------------------------------------------------------------------
c     Modified for use with the atomic energy network (aenet) package
c     2017-07-21 Nongnuch Artrith and Alexander Urban
c ----------------------------------------------------------------------
c
      subroutine extra1
c
      use aenettinker
c
      implicit none
c
c     Initialize aenet-Tinker interface
c     (Will only do something when called for first time)
      call aenet_tinker_init
c
c     Evaluate energy
      call aenet_eval_energy_and_forces
c
      return
      end
