module datastr
type :: node
    integer :: ndim
    integer, dimension(:), allocatable :: x
    type(node), pointer :: next => null()
end type node

type :: queue
    integer :: length
    type(node), pointer :: head => null()
    type(node), pointer :: tail => null()
end type queue

contains
subroutine node_new(pn, coords)
    type(node), pointer, intent(inout) :: pn
    integer, dimension(:), intent(in) :: coords
    
    if(associated(pn)) then
        nullify(pn)
    endif
    nullify(pn)
    allocate(pn)
    
    pn%ndim = size(coords)
    allocate(pn%x(size( coords )))
    pn%x = coords
end subroutine node_new

subroutine node_del(this)
    type(node), pointer, intent(inout) :: this
    deallocate(this%x)
    this%next => null()
end subroutine node_del

subroutine queue_init(this)
    type(queue), intent(inout) :: this
    this%length = 0
    this%head => null()
    this%tail => null()
end subroutine queue_init

subroutine queue_put(this, pnode)
    type(queue), intent(inout) :: this
    type(node), pointer, intent(inout) ::  pnode
    if(this%length .eq. 0) then
        this%head => pnode
        this%tail => pnode
        pnode%next => null()
    else
        this%tail%next => pnode
        pnode%next => null()
    endif
    pnode => null()
    this%length = this%length + 1
end subroutine queue_put

subroutine queue_get(this, pnode)
    type(queue), intent(inout) :: this
    type(node), pointer, intent(inout) ::  pnode
    if(this%length .eq. 0) then
        nullify(pnode)
        return
    elseif(this%length .eq. 1) then
        pnode => this%head
        this%head => null()
        this%tail => null()
        this%length = 0
    else
        pnode => this%head
        pnode%next => null()
        this%head => this%head%next
        this%length = this%length - 1
    endif
end subroutine queue_get

subroutine queue_del(this)
    type(queue), intent(inout) :: this
    type(node), pointer ::  pn
    if(this%length .eq. 0) then
        return
    else
        pn => this%head
        do while ( associated(pn) )
            pn => pn%next
            call node_del(this%head)
            this%head => pn
        enddo
    endif
end subroutine queue_del

subroutine queue_show(this)
    type(queue), intent(inout) :: this
    type(node), pointer ::  pn
    if(this%length .eq. 0) then
        return
    else
        print *, "queue"
        pn => this%head
        do while ( associated(pn) )
            pn => pn%next
            print *, this%head%x
            this%head => pn
        enddo
    endif
end subroutine queue_show

end module datastr

program test
use datastr
implicit none
    type(node), pointer :: pn1, pn2
    type(queue) :: pque
    integer, dimension(2) :: x = (/1,2/), y=(/2,3/)
    call node_new(pn1, x)
    pn2 => pn1
    call node_new(pn1, y)
    print *, pn1%x
    print *, pn2%x
    call queue_init(pque)
    call queue_put(pque,pn1)
    call queue_put(pque,pn2)
    x = x**2
    call node_new(pn1, x)
    call queue_put(pque,pn1)
    x = x**2
    call node_new(pn1, x)
    call queue_put(pque,pn1)
    call queue_show(pque)
end program test
