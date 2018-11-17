module datastr
type :: node
    integer :: ndim
    integer, dimension(:), allocatable :: x
    type(node), pointer :: next = null()
    contains
        procedure :: del => node_del
end type node

type :: queue
    integer :: length
    type(node), pointer :: head = null()
    type(node), pointer :: tail = null()
    contains
        procedure :: init => queue_init
        procedure :: put => queue_put
        procedure :: get => queue_get
        procedure :: del => queue_del 
end type queue

contains
subroutine node_cre(pn, coords)
    type(node), pointer, intent(inout) :: pn
    integer, dimension(:), intent(in) :: coords
    type(node), target :: node
    
    node%ndim = size(coords)
    allocate(node%x(size( coords )))
    node%x = coords
    pn => node
end subroutine node_cre

subroutine node_del(this)
    type(node), intent(inout) :: this
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
        inode%next => null()
    else
        this%tail%next => pnode
        inode%next => null()
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
        this%head => head%next
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
        do while (pn)
            pn => pn%next
            this%head%del()
            this%head => pn
        enddo
    endif
end subroutine queue_del

end module datastr

program test
use datastr
implicit none
    type(node), pointer :: pn
    type(queue) :: pque
    real(8), dimension(2) :: x = (/1,2/)
    call node_cre(pn, x)
    print *, pn%x
end program test
