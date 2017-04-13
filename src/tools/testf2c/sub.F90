module abcd
  implicit none

  integer :: NumInt111


end  module abcd

subroutine Sub_Fortran(NumInt,NumFloat,NumDouble)
  use abcd
  implicit none
  integer :: NumInt
  real :: NumFloat
  real(8) :: NumDouble

  NumInt111=10

  NumDouble=NumFloat**NumInt*NumInt111

  NumInt111=10

end subroutine

real(8) function Function_Fortran(NumDouble)
  use abcd

  implicit none
  real(8) :: NumDouble
  Function_Fortran = sqrt(NumDouble)*NumInt111

end function
