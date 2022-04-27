program main
  use globals
  implicit none
  call initialise()
  call logger%error("main","test")
end program main