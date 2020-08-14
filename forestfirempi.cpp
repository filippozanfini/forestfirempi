#include <mpi.h>
#include <allegro5/allegro.h>
#include <allegro5/allegro_primitives.h>
#include <stdio.h>
#include <math.h>

#define N 500

enum states
{
    EMPTY = 0,
    TREE = 1,
    FIRE = 2
};

float probRand()
{
    return (float)rand() / RAND_MAX;
}

bool withinGrid(int row, int col)
{
    if ((row < 0) || (col < 0))
        return false;
    if ((row >= N) || (col >= N))
        return false;

    return true;
}

bool getNeighbors(int x, int y, int **forest)
{
    for (int i = x - 1; i <= (x + 1); i++)
    {
        for (int j = y - 1; j <= (y + 1); j++)
        {
            if (!(i == x && j == y))
            {
                if (withinGrid(i, j) && forest[i][j] == FIRE)
                    return true;
            }
        }
    }

    return false;
}

void print(int **forest)
{
    al_clear_to_color(al_map_rgb(0, 0, 0));

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            switch (forest[i][j])
            {
            case EMPTY:
                al_draw_filled_rectangle(i * 5, j * 5, i * 5 + 5, j * 5 + 5, al_map_rgb(0, 0, 0));
                break;

            case TREE:
                al_draw_filled_rectangle(i * 5, j * 5, i * 5 + 5, j * 5 + 5, al_map_rgb(34, 139, 34));
                break;

            case FIRE:
                al_draw_filled_rectangle(i * 5, j * 5, i * 5 + 5, j * 5 + 5, al_map_rgb(255, 140, 0));
                break;
            }
        }
    }

    al_flip_display();
    al_rest(0.2);
}

void setRules(int **recvBuffer, int **forest, int i, int j, float p, float f)
{
    switch (recvBuffer[i][j])
    {
    case EMPTY:
        if (probRand() < p)
            forest[i][j] = TREE;
        break;

    case TREE:
        if (getNeighbors(i, j, recvBuffer) || probRand() < f)
            forest[i][j] = FIRE;
        break;

    case FIRE:
        forest[i][j] = EMPTY;
        break;

    default:
        forest[i][j] = recvBuffer[i][j];
        break;
    }
}

int main(int argc, char *argv[])
{

    int rank, size;
    MPI_Status status;
    MPI_Request request;
    MPI_Request sendLeftRqst;
    MPI_Request sendRightRqst;
    MPI_Request recvLeftRqst;
    MPI_Request recvRightRqst;

    int portion;

    float p = 0.02;   // Probability of new growth
    float f = 0.0002; // Probability of lightning

    int done = 1;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (N % (size - 1) != 0)
    {
        if (rank == 0)
            printf("Invalid number of processes\n");
        MPI_Finalize();
        return 0;
    }

    portion = N / (size - 1);

    // Buffer matrix
    int **rcvBuffer = new int *[portion + 2];
    for (int i = 0; i < portion + 2; i++)
    {
        rcvBuffer[i] = new int[N];
    }

    // Future generation matrix
    int **forest = new int *[N];
    for (int i = 0; i < N; i++)
    {
        forest[i] = new int[N];
    }

    if (rank == 0)
    {

        int **currentForest = new int *[N];
        for (int i = 0; i < N; i++)
        {
            currentForest[i] = new int[N];
        }

        // Allegro setup

        ALLEGRO_DISPLAY *display;
        const int altezzaSchermo = N;
        const int larghezzaSchermo = N;
        const int lato = 10;

        ALLEGRO_KEYBOARD_STATE key_state;

        al_init();
        display = al_create_display(larghezzaSchermo, altezzaSchermo);
        al_init_primitives_addon();

        if (!al_install_keyboard())
        {
            printf("Error with keyboard\n");
            MPI_Finalize();
            return -1;
        }

        // currentForest init
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                currentForest[i][j] = EMPTY;
            }
        }

        while (done != 0)
        {
            al_get_keyboard_state(&key_state);
            if (al_key_down(&key_state, ALLEGRO_KEY_Z))
            {
                done = 0;
            }
            for (int j = 1; j < size; j++)
                MPI_Send(&done, 1, MPI_INT, j, 99, MPI_COMM_WORLD);

            int k = 0;

            for (int j = 1; j < size; j++)
            {
                for (int i = 0; i < portion; i++)
                {
                    MPI_Isend(&(currentForest[i + k][0]), N, MPI_INT, j, 11, MPI_COMM_WORLD, &request);
                }
                k += portion;
            }

            // Print with Allegro
            print(currentForest);

            k = 0;
            for (int j = 1; j < size; j++)
            {
                for (int i = 0; i < portion; i++)
                {
                    MPI_Recv(&(currentForest[i + k][0]), N, MPI_INT, j, 75, MPI_COMM_WORLD, &status);
                }

                k += portion;
            }
        }

        al_destroy_display(display);

        for (int i = 0; i < N; i++)
        {
            delete[] currentForest[i];
        }
        delete[] currentForest;
    }

    else if (rank == 1)
    {
        while (done != 0)
        {

            MPI_Recv(&done, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);

            for (int i = 0; i < portion; i++)
            {
                MPI_Recv(&(rcvBuffer[i][0]), N, MPI_INT, 0, 11, MPI_COMM_WORLD, &status);
            }

            if ((size - 1) != 1)
            {
                MPI_Isend(&(rcvBuffer[portion - 1][0]), N, MPI_INT, rank + 1, 12, MPI_COMM_WORLD, &sendRightRqst);
                MPI_Irecv(&(rcvBuffer[portion][0]), N, MPI_INT, rank + 1, 13, MPI_COMM_WORLD, &recvRightRqst);
                MPI_Wait(&sendRightRqst, &status);
                MPI_Wait(&recvRightRqst, &status);
            }

            for (int i = 1; i < portion; i++)
                for (int j = 0; j < N; j++)
                    setRules(rcvBuffer, forest, i, j, p, f);

            for (int i = 0; i < portion; i++)
            {
                MPI_Send(&(forest[i][0]), N, MPI_INT, 0, 75, MPI_COMM_WORLD);
            }
        }
    }
    else if (rank == size - 1)
    {
        while (done != 0)
        {

            MPI_Recv(&done, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);

            for (int i = 1; i < portion + 1; i++)
            {
                MPI_Recv(&(rcvBuffer[i][0]), N, MPI_INT, 0, 11, MPI_COMM_WORLD, &status);
            }

            MPI_Isend(&(rcvBuffer[1][0]), N, MPI_INT, rank - 1, 13, MPI_COMM_WORLD, &sendLeftRqst);
            MPI_Irecv(&(rcvBuffer[0][0]), N, MPI_INT, rank - 1, 12, MPI_COMM_WORLD, &recvLeftRqst);

            MPI_Wait(&sendLeftRqst, &status);
            MPI_Wait(&recvLeftRqst, &status);

            for (int i = 1; i < portion + 1; i++)
                for (int j = 0; j < N; j++)
                    setRules(rcvBuffer, forest, i, j, p, f);

            for (int i = 1; i < portion + 1; i++)
            {
                MPI_Send(&(forest[i][0]), N, MPI_INT, 0, 75, MPI_COMM_WORLD);
            }
        }
    }
    else
    {
        while (done != 0)
        {

            MPI_Recv(&done, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);

            for (int i = 1; i < portion + 1; i++)
            {
                MPI_Recv(&(rcvBuffer[i][0]), N, MPI_INT, 0, 11, MPI_COMM_WORLD, &status);
            }

            MPI_Isend(&(rcvBuffer[1][0]), N, MPI_INT, rank - 1, 13, MPI_COMM_WORLD, &sendLeftRqst);
            MPI_Isend(&(rcvBuffer[portion][0]), N, MPI_INT, rank + 1, 12, MPI_COMM_WORLD, &sendRightRqst);
            MPI_Irecv(&(rcvBuffer[0][0]), N, MPI_INT, rank - 1, 12, MPI_COMM_WORLD, &recvLeftRqst);
            MPI_Irecv(&(rcvBuffer[portion + 1][0]), N, MPI_INT, rank + 1, 13, MPI_COMM_WORLD, &recvRightRqst);

            MPI_Wait(&sendLeftRqst, &status);
            MPI_Wait(&sendRightRqst, &status);
            MPI_Wait(&recvLeftRqst, &status);
            MPI_Wait(&recvRightRqst, &status);

            for (int i = 1; i < portion + 1; i++)
                for (int j = 0; j < N; j++)
                    setRules(rcvBuffer, forest, i, j, p, f);

            for (int i = 1; i < portion + 1; i++)
            {

                MPI_Send(&(forest[i][0]), N, MPI_INT, 0, 75, MPI_COMM_WORLD);
            }
        }
    }

    for (int i = 0; i < portion + 2; i++)
    {
        delete[] rcvBuffer[i];
    }
    delete[] rcvBuffer;

    for (int i = 0; i < N; i++)
    {
        delete[] forest[i];
    }
    delete[] forest;

    MPI_Finalize();
    return 0;
}