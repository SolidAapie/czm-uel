!
!      Return information about an element according to its type
!

      subroutine INT1_element_info_nope(element_type, nope)
!
!      Returns number of nodes per element
!
!      ARGUMENTS
!      element_type     string identifier
!
!      RETURNS
!      nope             nodes per element
      implicit none
      character*6 element_type
      integer nope

c$$$      if((element_type.eq.'C3D20R').or.(element_type.eq.'C3D20 ')) then
c$$$         nope=20
c$$$      elseif((element_type.eq.'C3D8R ').or.(element_type.eq.'C3D8  ')) 
c$$$     &        then
c$$$         nope=8
c$$$      elseif(element_type.eq.'INT1  ') then
c          nope=8
          nope=20
c$$$      else
c$$$         nope=10
c$$$      endif

      return
      end



      subroutine INT1_element_info_nopev(element_type, nopev)
!
!      Returns number of nodes at element vertices
!
!      ARGUMENTS
!      element_type     string identifier
!
!      RETURNS 
!      nopev             nodes per element vertex(?)
      implicit none
      character*6 element_type
      integer nopev

c$$$      if((element_type.eq.'C3D20R').or.(element_type.eq.'C3D20 ')) then
c$$$         nopev=8
c$$$      elseif((element_type.eq.'C3D8R ').or.(element_type.eq.'C3D8  ')) 
c$$$     &        then
c$$$         nopev=8
c$$$      elseif(element_type.eq.'INT1  ') 
c$$$     &        then
         nopev=8
c$$$      else
c$$$         nopev=4
c$$$      endif

      return
      end



      subroutine INT1_element_info_nopes(element_type, nopes)
!
!      Returns number of nodes on an element surface
!
!      ARGUMENTS
!      element_type     string identifier
!
!      RETURNS
!      nope             nodes per element
      implicit none
      character*6 element_type
      integer nopes

c$$$      if((element_type.eq.'C3D20R').or.(element_type.eq.'C3D20 ')) then
c$$$         nopes=8
c$$$      elseif((element_type.eq.'C3D8R ').or.(element_type.eq.'C3D8  ')) 
c$$$     &        then
c$$$         nopes=4
c$$$c$$$      elseif(element_type.eq.'INT1  ') 
c$$$c$$$     &        then
c$$$c$$$         nopes=4
c$$$      else
         nopes=6
c$$$      endif

      return
      end



      subroutine INT1_element_info_mint2d(element_type, mint2d)
!
!      Returns number of integration points (or, more generally,
!      information points) in 2d
!
!      ARGUMENTS
!      element_type     string identifier
!
!      RETURNS
!      mint2d           number of integration points in 2d
      implicit none
      character*6 element_type
      integer mint2d

c$$$      if(element_type.eq.'C3D8R ') then
c$$$         mint2d=1
c$$$      elseif((element_type.eq.'C3D8  ').or.(element_type.eq.'C3D20R'))
c$$$     &        then
c$$$         mint2d=4
c$$$      elseif(element_type.eq.'C3D20 ') then
c$$$         mint2d=9
c$$$c$$$      elseif(element_type.eq.'INT1  ') then
c$$$c$$$         mint2d=0
c$$$      else
         mint2d=3
c$$$      endif

      return
      end



      subroutine INT1_element_info_mint3d(element_type, mint3d)
!
!      Returns number of integration points (or, more generally,
!      information points) in 3d
!
!      ARGUMENTS
!      element_type     string identifier
!
!      RETURNS
!      mint3d           number of integration points in 2d
      implicit none
      character*6 element_type
      integer mint3d

c$$$      if(element_type.eq.'C3D8R ') then
c$$$         mint3d=1
c$$$      elseif((element_type.eq.'C3D8  ').or.(element_type.eq.'C3D20R'))
c$$$     &        then
c$$$         mint3d=8
c$$$      elseif(element_type.eq.'C3D20 ') then
c$$$         mint3d=27
c$$$      elseif(element_type.eq.'INT1  ') then
         mint3d=9
c$$$      else
c$$$         mint3d=4
c$$$      endif

      return
      end



