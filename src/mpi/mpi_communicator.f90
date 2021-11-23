module mpi_communicator
    use mpi_f08
    implicit none

    type(MPI_Comm), parameter :: comm = MPI_COMM_WORLD
    integer                   :: mpi_err
    integer                   :: mpi_rank
    integer                   :: mpi_size

    contains

        subroutine mpi_comm_initialise
            call MPI_Init(mpi_err)
            call MPI_Comm_size(comm, mpi_size, mpi_err)
            call MPI_Comm_rank(comm, mpi_rank, mpi_err)
            call MPI_Comm_size(comm, mpi_size, mpi_err)
        end subroutine mpi_comm_initialise

        subroutine mpi_comm_finalise
            logical :: flag
            call MPI_Initialized(flag, mpi_err)

            if (flag) then
                call MPI_Finalize(mpi_err)
            endif

            call MPI_Finalized(flag, mpi_err)

            if (.not. flag) then
                call MPI_Abort(comm, -1, mpi_err)
            endif
        end subroutine mpi_comm_finalise
end module mpi_communicator
