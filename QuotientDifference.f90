Program QuotientDifference

!1- Ingresa el grado del polinomio
!2- Ingresa los coeficientes de forma ordenada
!3- Comienza QD
    !i-Ejecuta primeras iteraciones

IMPLICIT NONE

integer:: Grado, i ,Metodo
real(8), allocatable, dimension(:) :: VecCoeficientes,VecRaices
real(8) tol

! COMIENZA EL PROGRAMA

CALL LeeCoeficientes(VecCoeficientes,Grado)
!----------ERROR--------!
write(*,*)"Ingrese la cota maxima de error" 
read(*,*)tol
!------Iteraciones-----!
write(*,*)"Ingrese maxima cantidad de iteraciones" 
read(*,*)maxiter

Metodo=0
CALL PreProcesamiento(VecCoeficientes,Grado,Metodo,Diferencia)
do while (CoeficienteCero(VecCoeficientes,Grado) .AND. Metodo<4) 

	Metodo=Metodo+1
	CALL PreProcesamiento(VecCoeficientes,Grado,Metodo,Diferencia) !hay 3 metodos para corregir 0s, haciendo diferentes reempazos s=x+d s+cx s=1/x
	!ademas puede bajar de grado si el ultimo termino es 0
	
enddo

if (CoeficienteCero(VecCoeficientes,Grado)) then !si el metodo es 4 enotnces no fue posible alterar el polinomio

  write(*,*)"Imposible obtener raices del polinomio : No fue posible modificarlo para que sus coeficientes fueran distintos de 0")
  
 else
	select case (grado)
	case (1) !grado=0
		write(*,*)"Error grado=0"
	case (2)
		write(*,*)"Raiz: "VecCoeficientes(1)/VecCoeficientes(2)
	case (3)
		call Resolvente(VecCoeficientes)!Falta desarrollar
	case default
		call QuotientDifferenceSub(VecCoeficientes,VecRaices,grado,tol,maxiter)
		call CorrigeSolucion(VecRaices,grado,Metodo,diferencia)
	end select
endif




CONTAINS

!----------------------------------!
subroutine PreProcesamiento(VecCoeficientes,Grado,Metodo,Diferencia)
real(8) ,intent(inout) :: VecCoeficientes(:)
real(8) intent(out) :: Diferencia
integer, intent(inout) ::Metodo,Grado


select case(Metodo)
	case (0) !Caso 0 en ultimo lugar
		do while (VecCoeficientes(grado)==0)
			grado=grado-1
		enddo
	case (1)
		call Metodo1(VecCoeficientes,Grado)
	case (2)
	case (3)
end select




endsubroutine

subroutine Metodo1(VecCoeficientes,Grado)
real(8) ,intent(inout) :: VecCoeficientes(:)
integer, intent(inout) :: Grado
integer i,j

do i=1,grado
	do j=1,grado
		VecCoeficientes(i)=Combinatoria(j,grado)*d*VecCoeficientes!esta mal
	enddo
enddo

endsubroutine


function Combinatoria(m,n)
integer intent(in) :: m,n
real(8) Combinatoria

Combinatoria=Factorial(m)/(Factorial(n)*Factorial(m-n)) 

endfunction

function Factorial(x)
integer, intent(in) :: x
integer Factorial

if x==0 then
	Factorial=1
else
	Factorial=1
	do i=1,x
		Factorial=Factorial*i
		
	enddo
endfunction
!----------------------------------!

function CoeficienteCero(VecCoeficientes,Grado)
logical CoeficienteCero
real(8) VecCoeficientes(:)
integer Grado,i

do while ((VecCoeficientes(i)/=0) .AND. (i<Grado))
	i=i+1
enddo

CoeficienteCero=((Grado==i).AND.(VecCoeficientes(grado)/=0))

endfunction

!----------------------------------!
SUBROUTINE LeeCoeficientes(VecCoeficientes,Grado)
integer i,Grado
real(8), allocatable :: VecCoeficientes(:)

open(2,FILE="CoeficientesPolinomio.txt", action="read")
read(2,*)grado
grado=grado+1
allocate(VecCoeficientes(grado))

do i=1,grado 
	read(2,*)VecCoeficientes(i)
enddo



END SUBROUTINE
!----------------------------------!

subroutine QuotientDifferenceSub(VecCoeficientes,VecRaices,tol,maxiter)
!parametros 
real(8) intent(inout)::VecCoeficientes(:),VecRaices(:)
integer ,intent(in):: grado,maxiter
real(8) ,intent(in)::tol
!var locales
integer i
real(8) error
real(8), allocatable, dimension(:) :: VecError
!VECE va de 1 a grado
!VEcRaiz va de 1 a grado
!VecCoeficientes va de 1 a grado
allocate(VecError(grado))

   
!COMIENZA QD
CALL CondicionesIniciales(VecCoeficientes,VecRaices,VecError,grado) 

iter=1
error=2*tol !por que?

do while (error>tol .AND. maxiter>iter)

    do i=1, grado                  !Calculo la fila de entera en cada iteración y el grado +1 porque el vector comiena en 1 , el el coef de maximo grado esta een grado +1
       VecRaices(i)=VecError(i)-VecError(i-1)+VecRaices(i)     
    enddo
    
    Do i=2, grado-1            					    !Como e(0)=e(n)=0 siempre, uso esos límites en el ciclo
       VecError(i)=VecError(i)*(VecRaices(i+1)/VecRaices(i))   ! Hay un error de tipeo en el pdf que figura q(i-1) pero en la tabla utiliza q(i+1), sino no da
    end do
   

     iter=iter+1
     error=maxval(abs(e))
end do

end subroutine

!----------------------------------!

subroutine CondicionesIniciales(VecCoeficientes,VecRaices,VecError,grado)
real(8) intent(inout)::VecCoeficientes(:),VecRaices(:)
integer ,intent(in):: grado,maxiter
real(8) intent(inout) :: VecError(:)

VecRaices=0
VecRaices(1)=-(VecCoeficientes(n-1))/(VecCoeficientes(n))

do i=2, n-1

      VecError(i)=VecCoeficientes(grado-i-1)/VecCoeficientes(grado-i)
      
end do

VecError(0)=0
VecError(grado)=0

endsubroutine
!----------------------------------!
End program
