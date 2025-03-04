'''Helper functions'''

def is_on_board(y, x, board):
    #Check if the position (y, x) is within the bounds of the board.
    return 0 <= y < len(board) and 0 <= x < len(board[0])

def check_sequence(board, y, x, d_y, d_x):
    colour = board[y][x]

    if colour == " ":  # No sequence can start with an empty space
        return None

    # Traverse the direction to check if a sequence of the same color exists
    for step in range(1, 5):
        ny, nx = y + step * d_y, x + step * d_x
        if not is_on_board(ny, nx, board) or board[ny][nx] != colour:
            return None  # Sequence is broken

    return colour  # Valid sequence of exactly '5'



def is_empty(board):

    #Get board dimensions
    board_size = len(board)

    for r in range(board_size):
        for c in range(board_size):
            if board[r][c] != " ":
                return False
    return True


def is_bounded(board, y_end, x_end, length, d_y, d_x):
    # Get starting point of the sequence
    y_start = y_end - d_y * (length - 1)
    x_start = x_end - d_x * (length - 1)

    # Check the space before the sequence
    y_before = y_start - d_y
    x_before = x_start - d_x
    if is_on_board(y_before, x_before, board):
        before = board[y_before][x_before]
    else:
        before = "BLOCKED"  # Out of board bounds

    # Check the space after the sequence
    y_after = y_end + d_y
    x_after = x_end + d_x
    if is_on_board(y_after, x_after, board):
        after = board[y_after][x_after]
    else:
        after = "BLOCKED"  # Out of board bounds

    # Determine the status of the sequence
    if before == " " and after == " ":
        '''print(f"Before: {before} After: {after}")
        print("open")'''
        return "OPEN"

    if (before == " " and after == "BLOCKED") or (before == " " and after in ["b", "w"]):
        '''print(f"Before: {before} After: {after}")
        print("semiopen")'''
        return "SEMIOPEN"

    if (after == " " and before == "BLOCKED") or (after == " " and before in ["b", "w"]):
        '''print(f"Before: {before} After: {after}")
        print("semiopen")'''
        return "SEMIOPEN"

    '''print(f"Before: {before} After: {after}")
    print("closed")'''
    return "CLOSED"


def detect_row(board, col, y_start, x_start, length, d_y, d_x):
    # Dimensions of the board
    max_y, max_x = len(board), len(board[0])

    # Counters for open and semi-open sequences
    open_seq_count = 0
    semi_open_seq_count = 0

    # Starting coordinates for search
    y, x = y_start, x_start

    # Iterate along the specified direction until the board boundary is reached
    while 0 <= y < max_y and 0 <= x < max_x:

        # Check if the current square contains the target color
        if board[y][x] == col:
            seq_len = 1
            is_longer = False  # Track if the sequence is longer than needed

            # Traverse in the given direction to calculate sequence length
            for step in range(1, length + 1):  # Check up to (length + 1) steps
                ny, nx = y + step * d_y, x + step * d_x

                # Check if within bounds and the color matches
                if 0 <= ny < max_y and 0 <= nx < max_x and board[ny][nx] == col:
                    seq_len += 1
                    if seq_len > length:  # If sequence is too long, mark it
                        is_longer = True
                        break
                else:
                    break

            # Only process if the sequence matches the exact target length
            if seq_len == length and not is_longer:
                # Coordinates of the last stone in the sequence
                ny, nx = y + (length - 1) * d_y, x + (length - 1) * d_x

                # Check if the sequence is open, semi-open, or closed
                status = is_bounded(board, ny, nx, length, d_y, d_x)

                # Update counts based on the sequence status
                if status == "OPEN":
                    open_seq_count += 1
                elif status == "SEMIOPEN":
                    semi_open_seq_count += 1

            # Skip to the end of the current sequence to avoid double counting
            y += seq_len * d_y
            x += seq_len * d_x
        else:
            # Move to the next square in the direction
            y += d_y
            x += d_x

    return open_seq_count, semi_open_seq_count



def detect_rows(board, col, length):
    open_seq_count, semi_open_seq_count = 0, 0


    # Check horizontal and diagonal sequences from the left edge
    for left_step in range(len(board)):
        # Horizontal sequences
        temp_open, temp_semi_open = detect_row(board, col, left_step, 0, length, 0, 1)
        open_seq_count += temp_open
        semi_open_seq_count += temp_semi_open

        # Upper-left to lower-right diagonals
        temp_open, temp_semi_open = detect_row(board, col, left_step, 0, length, 1, 1)
        open_seq_count += temp_open
        semi_open_seq_count += temp_semi_open

    # Check vertical sequences from the bottom edge
    # Check diagonals from the top edge (avoiding top-left corner)
    for top_step in range(1, len(board[0])):

        # Vertical Sequences
        temp_open, temp_semi_open = detect_row(board, col, 0, top_step, length, 1, 0)
        open_seq_count += temp_open
        semi_open_seq_count += temp_semi_open

        # Upper-left to lower-right diagonals
        temp_open, temp_semi_open = detect_row(board, col, 0, top_step, length, 1, 1)
        open_seq_count += temp_open
        semi_open_seq_count += temp_semi_open

        # Upper-right to lower-left diagonals
        temp_open, temp_semi_open = detect_row(board, col, 0, top_step, length, 1, -1)
        open_seq_count += temp_open
        semi_open_seq_count += temp_semi_open

    # Check diagonals from the right edge
    for right_step in range(1, len(board)):
        # Upper-right to lower-left diagonals
        temp_open, temp_semi_open = detect_row(board, col, right_step, 7, length, 1, -1)
        open_seq_count += temp_open
        semi_open_seq_count += temp_semi_open

    return open_seq_count, semi_open_seq_count


def search_max(board):

    max_score = -100001
    move_y, move_x = 0, 0

    # Place a black piece on every unoccupied tile on the board
    for y in range(len(board)):
        for x in range(len(board[0])):
            if board[y][x] == " ":
                board[y][x] = "b"

                # Check the score with this new board state
                if score(board) > max_score:
                    move_y, move_x = y, x
                    max_score = score(board)

                # Return board to original state
                board[y][x] = " "


    return move_y, move_x

''' Detects and returns the weight of each move '''
def score(board):
    MAX_SCORE = 100000

    open_b = {}
    semi_open_b = {}
    open_w = {}
    semi_open_w = {}

    for i in range(2, 6):
        open_b[i], semi_open_b[i] = detect_rows(board, "b", i)
        open_w[i], semi_open_w[i] = detect_rows(board, "w", i)


    if open_b[5] >= 1 or semi_open_b[5] >= 1:
        return MAX_SCORE

    elif open_w[5] >= 1 or semi_open_w[5] >= 1:
        return -MAX_SCORE

    return (-10000 * (open_w[4] + semi_open_w[4])+
            500  * open_b[4]                     +
            50   * semi_open_b[4]                +
            -100  * open_w[3]                    +
            -30   * semi_open_w[3]               +
            50   * open_b[3]                     +
            10   * semi_open_b[3]                +
            open_b[2] + semi_open_b[2] - open_w[2] - semi_open_w[2])


def is_win(board):
    max_y, max_x = len(board), len(board[0])
    available_spot = False

    for y in range(max_y):
        for x in range(max_x):
            if(
            check_sequence(board, y, x, 1, 0) or
            check_sequence(board, y, x, 0, 1) or
            check_sequence(board, y, x, 1, 1) or
            check_sequence(board, y, x, 1, -1)
            ) == "b":
                return "Black won"

            if(
            check_sequence(board, y, x, 1, 0) or
            check_sequence(board, y, x, 0, 1) or
            check_sequence(board, y, x, 1, 1) or
            check_sequence(board, y, x, 1, -1)
            ) == "w":
                return "White won"

            if board[y][x] == " ":
                available_spot = True

    if not available_spot:
        return "Draw"
    else:
        return "Continue playing"




def print_board(board):

    s = "*"
    for i in range(len(board[0])-1):
        s += str(i%10) + "|"
    s += str((len(board[0])-1)%10)
    s += "*\n"

    for i in range(len(board)):
        s += str(i%10)
        for j in range(len(board[0])-1):
            s += str(board[i][j]) + "|"
        s += str(board[i][len(board[0])-1])

        s += "*\n"
    s += (len(board[0])*2 + 1)*"*"

    print(s)


def make_empty_board(sz):
    board = []
    for i in range(sz):
        board.append([" "]*sz)
    return board



def analysis(board):
    for c, full_name in [["b", "Black"], ["w", "White"]]:
        print("%s stones" % (full_name))
        for i in range(2, 6):
            open, semi_open = detect_rows(board, c, i);
            print("Open rows of length %d: %d" % (i, open))
            print("Semi-open rows of length %d: %d" % (i, semi_open))

def play_gomoku(board_size):
    board = make_empty_board(board_size)
    board_height = len(board)
    board_width = len(board[0])

    while True:
        print_board(board)
        if is_empty(board):
            move_y = board_height // 2
            move_x = board_width // 2
        else:
            move_y, move_x = search_max(board)

        print("Computer move: (%d, %d)" % (move_y, move_x))
        board[move_y][move_x] = "b"
        print_board(board)
        analysis(board)

        game_res = is_win(board)
        if game_res in ["White won", "Black won", "Draw"]:
            return game_res


        print("Your move:")
        move_y = int(input("y coord: "))
        move_x = int(input("x coord: "))
        board[move_y][move_x] = "w"
        print_board(board)
        analysis(board)

        game_res = is_win(board)
        if game_res in ["White won", "Black won", "Draw"]:
            return game_res


def put_seq_on_board(board, y, x, d_y, d_x, length, col):
    for i in range(length):
        board[y][x] = col
        y += d_y
        x += d_x


def test_is_empty():
    board  = make_empty_board(8)
    if is_empty(board):
        print("TEST CASE for is_empty PASSED")
    else:
        print("TEST CASE for is_empty FAILED")

def test_is_bounded():
    board = make_empty_board(8)
    x = 5; y = 1; d_x = 0; d_y = 1; length = 3
    put_seq_on_board(board, y, x, d_y, d_x, length, "w")
    print_board(board)

    y_end = 3
    x_end = 5

    if is_bounded(board, y_end, x_end, length, d_y, d_x) == 'OPEN':
        print("TEST CASE for is_bounded PASSED")
    else:
        print("TEST CASE for is_bounded FAILED")


def test_detect_row():
    board = make_empty_board(8)
    x = 5; y = 1; d_x = 0; d_y = 1; length = 3
    put_seq_on_board(board, y, x, d_y, d_x, length, "w")
    print_board(board)
    if detect_row(board, "w", 0,x,length,d_y,d_x) == (1,0):
        print("TEST CASE for detect_row PASSED")
    else:
        print("TEST CASE for detect_row FAILED")

def test_detect_rows():
    board = make_empty_board(8)
    x = 5; y = 1; d_x = 0; d_y = 1; length = 3; col = 'w'
    put_seq_on_board(board, y, x, d_y, d_x, length, "w")
    print_board(board)
    if detect_rows(board, col,length) == (1,0):
        print("TEST CASE for detect_rows PASSED")
    else:
        print("TEST CASE for detect_rows FAILED")

def test_search_max():
    board = make_empty_board(8)
    x = 5; y = 0; d_x = 0; d_y = 1; length = 4; col = 'w'
    put_seq_on_board(board, y, x, d_y, d_x, length, col)
    x = 6; y = 0; d_x = 0; d_y = 1; length = 4; col = 'b'
    put_seq_on_board(board, y, x, d_y, d_x, length, col)
    print_board(board)
    if search_max(board) == (4,6):
        print("TEST CASE for search_max PASSED")
    else:
        print("TEST CASE for search_max FAILED")

def easy_testset_for_main_functions():
    test_is_empty()
    test_is_bounded()
    test_detect_row()
    test_detect_rows()
    test_search_max()

def some_tests():
    board = make_empty_board(8)

    board[0][5] = "w"
    board[0][6] = "b"
    y = 5; x = 2; d_x = 0; d_y = 1; length = 3
    put_seq_on_board(board, y, x, d_y, d_x, length, "w")
    print_board(board)
    analysis(board)

    # Expected output:
    #       *0|1|2|3|4|5|6|7*
    #       0 | | | | |w|b| *
    #       1 | | | | | | | *
    #       2 | | | | | | | *
    #       3 | | | | | | | *
    #       4 | | | | | | | *
    #       5 | |w| | | | | *
    #       6 | |w| | | | | *
    #       7 | |w| | | | | *
    #       *****************
    #       Black stones:
    #       Open rows of length 2: 0
    #       Semi-open rows of length 2: 0
    #       Open rows of length 3: 0
    #       Semi-open rows of length 3: 0
    #       Open rows of length 4: 0
    #       Semi-open rows of length 4: 0
    #       Open rows of length 5: 0
    #       Semi-open rows of length 5: 0
    #       White stones:
    #       Open rows of length 2: 0
    #       Semi-open rows of length 2: 0
    #       Open rows of length 3: 0
    #       Semi-open rows of length 3: 1
    #       Open rows of length 4: 0
    #       Semi-open rows of length 4: 0
    #       Open rows of length 5: 0
    #       Semi-open rows of length 5: 0

    y = 3; x = 5; d_x = -1; d_y = 1; length = 2

    put_seq_on_board(board, y, x, d_y, d_x, length, "b")
    print_board(board)
    analysis(board)

    # Expected output:
    #        *0|1|2|3|4|5|6|7*
    #        0 | | | | |w|b| *
    #        1 | | | | | | | *
    #        2 | | | | | | | *
    #        3 | | | | |b| | *
    #        4 | | | |b| | | *
    #        5 | |w| | | | | *
    #        6 | |w| | | | | *
    #        7 | |w| | | | | *
    #        *****************
    #
    #         Black stones:
    #         Open rows of length 2: 1
    #         Semi-open rows of length 2: 0
    #         Open rows of length 3: 0
    #         Semi-open rows of length 3: 0
    #         Open rows of length 4: 0
    #         Semi-open rows of length 4: 0
    #         Open rows of length 5: 0
    #         Semi-open rows of length 5: 0
    #         White stones:
    #         Open rows of length 2: 0
    #         Semi-open rows of length 2: 0
    #         Open rows of length 3: 0
    #         Semi-open rows of length 3: 1
    #         Open rows of length 4: 0
    #         Semi-open rows of length 4: 0
    #         Open rows of length 5: 0
    #         Semi-open rows of length 5: 0
    #

    y = 5; x = 3; d_x = -1; d_y = 1; length = 1
    put_seq_on_board(board, y, x, d_y, d_x, length, "b");
    print_board(board);
    analysis(board);

    #        Expected output:
    #           *0|1|2|3|4|5|6|7*
    #           0 | | | | |w|b| *
    #           1 | | | | | | | *
    #           2 | | | | | | | *
    #           3 | | | | |b| | *
    #           4 | | | |b| | | *
    #           5 | |w|b| | | | *
    #           6 | |w| | | | | *
    #           7 | |w| | | | | *
    #           *****************
    #
    #
    #        Black stones:
    #        Open rows of length 2: 0
    #        Semi-open rows of length 2: 0
    #        Open rows of length 3: 0
    #        Semi-open rows of length 3: 1
    #        Open rows of length 4: 0
    #        Semi-open rows of length 4: 0
    #        Open rows of length 5: 0
    #        Semi-open rows of length 5: 0
    #        White stones:
    #        Open rows of length 2: 0
    #        Semi-open rows of length 2: 0
    #        Open rows of length 3: 0
    #        Semi-open rows of length 3: 1
    #        Open rows of length 4: 0
    #        Semi-open rows of length 4: 0
    #        Open rows of length 5: 0
    #        Semi-open rows of length 5: 0



if __name__ == '__main__':
    # To play with an 8x8 board
    play_gomoku(8)
