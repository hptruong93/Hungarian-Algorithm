'''
    Hungarian algorithm implementation. Details about the algorithm: https://en.wikipedia.org/wiki/Hungarian_algorithm
'''
import numpy
import itertools
import random

def debug(message):
    print message

def info(message):
    print message

_INFINITY = 99999999
_EPSILON = 0.0000001

_equal_zero = lambda x : abs(x) < _EPSILON

def _assertion(x, defensive = True):
    if defensive:
        assert x

def find_optimal_assignment(matrix, expect_optimal = True):
    '''
        Return a tuple of tuples (x,y) indicating the optimal assignment.
        Input:
            SQUARE matrix of size n x n
            expect_optimal: if this is True, the function asserts at the end that indeed n results are found
    '''
    #Dimension of the square matrix
    n = matrix.shape[0]
    output = []
    min_value = {
        'min_value' : _INFINITY,
        'min_coordinate' : None
    }

    # debug("Finding optimal solution...")
    # debug(matrix)

    zero_counts = tuple(count.tolist()[0] for count in _equal_zero(matrix).nonzero())

    #List of all coordinates of zeroes in the matrix
    zero_coordinates = list(zip(*zero_counts))

    def _pick(axis, value, index):
        '''
            Pick a coordinate as part of the solution
            axis is either 0 (x) or 1 (y)
            value is the index of the coordinate on the axis
            index is the index of the other axis on the other axis
        '''
        coordinate = (index, value) if axis == 1 else (value, index)
        # debug('Picking %s' % str(coordinate))
        _assertion(coordinate in zero_coordinates)

        min_value['min_coordinate'] = None
        min_value['min_value'] = _INFINITY

        #Ignore zeroes on the picked row and columns
        filtered_coordinates = [value for value in zero_coordinates if value[0] != coordinate[0] and value[1] != coordinate[1]]
        zero_coordinates[:] = [] #Empty the list
        zero_coordinates.extend(filtered_coordinates) #Regenerate the list of zeros in the matrix
        
        # debug('New zero list is %s' % zero_coordinates)
        output.append(coordinate)

    def _find_by_axis(axis, count_by_axis):
        '''
            Find the index of the row/col which has the least amount of zeros in the matrix.
            Will pick the entry as part of solution if the row/col has only one zero.
            Otherwise update the min_value variable
            Input
                axis is either 0 (x) or 1 (y)
                An iterator where each element is a tuple of 2 elements:
                    the index of the zero on the axis, and
                    an iterator containing index of the zeros on the other axis
        '''
        # debug("Finding with axis %s" % ('x' if axis == 0 else 'y'))
        for value, iterator in count_by_axis:
            zero_indices = list(iterator)
            # debug('For this axis %s, value is %s and zero indices are %s' % (axis, value, zero_indices))
            # debug('Index %s has %s indices' % (value, len(zero_indices)))

            if len(zero_indices) == 1:
                _pick(axis, value, zero_indices[0][1 - axis])
                return True
            else:
                # debug("This has %s zeros on it and min val is %s" % (len(zero_indices), min_value['min_value']))
                if min_value['min_value'] > len(zero_indices):
                    # debug("Found min %s with value %s on axis %s" % (len(zero_indices), value, 'x' if axis == 0 else 'y'))
                    min_value['min_value'] = len(zero_indices)
                    min_value['min_coordinate'] = (axis, value, zero_indices[0][1 - axis])

        return False

    def _find_unique_zero():
        '''
            Among the zeros coordinates, attempt to find the one that is unique on its row/col (i.e. it is the only zero on the row/col)
            If such zero exists, add it to the optimal solution
            If not, add the zero whose row/col has the smallest amount of zeros to the optimal solution
        '''
        for axis in range(len(zero_counts)):
            key_function = lambda pair : pair[axis]
            count_by_axis = itertools.groupby(sorted(zero_coordinates, key = key_function), key = key_function)
            if _find_by_axis(axis, count_by_axis):        
                return True

        return False

    while True:
        if not _find_unique_zero():
            _pick(*min_value['min_coordinate'])
            
        if len(zero_coordinates) == 0:
            break        

    # debug("Found %s" % output)
    _assertion(len(output) == n, expect_optimal)
    return output
    

def find_zero_covering_lines(matrix):
    '''
        Find the minimum number of horizontal and vertical lines that can cover all zeroes in the matrix
        Input:
            a square matrix with either
                at least one zero in each row or
                at least one zero in each column
        Return a tuple of two elements: (x_lines, y_lines) where
            x_lines is an iterable with indices of all covered rows
            y_lines is an iterable with indices of all covered columns
    '''
    _assertion(matrix.shape[0] == matrix.shape[1])
    n = matrix.shape[0]
    output = (set(), set())

    # debug("Finding zero covering lines...")
    # debug(matrix)

    zero_counts = tuple(count.tolist()[0] for count in _equal_zero(matrix).nonzero())
    
    #List of all coordinates of zeroes in the matrix
    zero_coordinates = sorted(list(zip(*zero_counts)))
    # debug("There are %s zeros" % len(zero_coordinates))

    using_matrix = numpy.matrix(matrix)

    marked_rows, marked_cols = set(), set()
    assigned_rows, assigned_cols, assignments = set(), set(), set()

    def _mark(x, y):
        # debug("Mark %s, %s" % (x, y))
        assigned_rows.add(x)
        assigned_cols.add(y)
        assignments.add((x, y))

    for x, y in find_optimal_assignment(using_matrix, expect_optimal = False):
        _mark(x, y)
    
    #Mark all rows have no assignment
    marked_rows = set(range(n)) - assigned_rows
    newly_marked_rows = set(marked_rows)

    while True:
        #Mark all columns (that haven't been marked) having zero in newly marked rows
        newly_marked_cols = set()
        for row in newly_marked_rows:
            new_cols = [col for col in xrange(n) if _equal_zero(matrix[row, col]) and col not in marked_cols]
            for col in new_cols:
                newly_marked_cols.add(col)
                marked_cols.add(col)

        if len(newly_marked_cols) == 0:
            break

        #Mark all rows (that haven't been marked) having assignment in newly marked cols
        newly_marked_rows = set()
        for col in newly_marked_cols:
            new_rows = [row for row in xrange(n) if (row, col) in assignments and row not in marked_rows]

            for row in new_rows:
                newly_marked_rows.add(row)
                marked_rows.add(row)

        if len(newly_marked_rows) == 0:
            break

    #Draw lines through all UNMARKED rows and marked columns
    output = (tuple(set(range(n)) - marked_rows), tuple(marked_cols))
    _assertion(sum(map(len, output)) <= len(zero_coordinates)) #Number of lines drawn is at most the number of zeros present

    # debug("Found %s. That's %s lines" % (str(output), sum(map(len, output))))
    return output

def smallest_uncovered_entry(matrix, covering_lines):
    '''
        Retrieve the smallest entry in the matrix that are not covered by the covering_lines
        Input:
            square matrix input
            a tuple (x_lines, y_lines) where
                x_lines is an iterable with indices of all covered rows
                y_lines is an iterable with indices of all covered columns

    '''
    # debug("Finding smallest uncovered entry for")
    # debug(matrix)
    # debug(covering_lines)

    x_lines, y_lines = covering_lines
    truncated_matrix = numpy.delete(matrix, x_lines, 0)
    truncated_matrix = numpy.delete(truncated_matrix, y_lines, 1)

    output = truncated_matrix.min()
    return output

def hungarian_algorithm_square(input_square_matrix):
    '''
        Apply hungarian algorithm to the input matrix
        Input:
            A SQUARE matrix
        Output:
            a tuple of tuples (x,y) indicating the optimal assignment.
    '''
    _assertion(input_square_matrix.shape[0] == input_square_matrix.shape[1])
    n = input_square_matrix.shape[0]
    using_matrix = numpy.matrix(input_square_matrix)

    # debug("Input matrix")
    # debug(input_square_matrix)

    using_matrix = using_matrix - using_matrix.min(1) #For each row, subtract that row by its min
    using_matrix = using_matrix - using_matrix.min(0) #For each column, subtract that col by its min

    # debug("After cleaning")
    # debug(using_matrix)

    while True:
        x_lines, y_lines = find_zero_covering_lines(using_matrix)
        _assertion(x_lines is not None and y_lines is not None)

        line_count = len(x_lines) + len(y_lines)
        _assertion(line_count <= n)
        if line_count == n:
            return find_optimal_assignment(using_matrix)
        else:
            smallest = smallest_uncovered_entry(using_matrix, (x_lines, y_lines))

        #Subtract the smallest value from all uncovered rows
        for row in [r for r in xrange(n) if r not in x_lines]:
            using_matrix[row] -= smallest

        #Add the smallest value to all covered columns
        for col in [c for c in xrange(n) if c in y_lines]:
            using_matrix[:, col] += smallest

def hungarian_algorithm_general(input_matrix):
    '''
        Apply hungarian algorithm to the input matrix
        Input:
            A rectangular matrix (can be square or non-square)
        Output:
            a tuple of tuples (x,y) indicating the optimal assignment.
    '''
    m = input_matrix.shape[0]
    n = input_matrix.shape[1]

    using_matrix = None
    if m > n: #Pad with 0s columns
        null_matrix = numpy.matrix(numpy.zeros((m, m - n)))
        using_matrix = numpy.append(input_matrix, null_matrix, 1)
    elif n > m: #Pad with 0s rows
        null_matrix = numpy.matrix(numpy.zeros((n - m, n)))
        using_matrix = numpy.append(input_matrix, null_matrix, 0)
    else: #Already square, doing nothing
        using_matrix = input_matrix

    return tuple((x, y) for x, y in hungarian_algorithm_square(using_matrix) if x < m and y < n)

def hungarian_algorithm_interface(iterable_of_iterable):
    '''
        Interface to this module. Apply hungarian algorithm to the input
        Input:
            An iterable of iterable. Has to be rectangular (i.e. each row has the same number of elements)
        Output:
            a tuple of tuples (x,y) indicating the optimal assignment.
    '''
    if not iterable_of_iterable:
        return None

    _assertion(len(set(map(len, iterable_of_iterable))) == 1) #assert rectangular iterable
    return hungarian_algorithm_general(numpy.matrix(iterable_of_iterable))

def termination_test(iterable_of_iterable):
    # debug("*" * 50)
    # debug(iterable_of_iterable)
    debug(hungarian_algorithm_interface(iterable_of_iterable))

def brute_force_solution(iterable_of_iterable):
    n = len(iterable_of_iterable)
    def _choices_generator(n):
        for permutation in itertools.permutations(xrange(n), n):
            yield zip(xrange(n), permutation)

    def _calculate_cost(choice):
        cost = sum(iterable_of_iterable[x][y] for x, y in choice)
        if cost < _calculate_cost.best_cost:
            _calculate_cost.best_choice = choice
            _calculate_cost.best_cost = cost

    _calculate_cost.best_choice = None
    _calculate_cost.best_cost = 99999999

    map(_calculate_cost, _choices_generator(n))
    return (_calculate_cost.best_choice, _calculate_cost.best_cost)

if __name__ == "__main__":
    pass
    # info("Start testing...")
    # xs = [[[1,2],[4,3]],
    #     [[40,0,5,0], [0,25,0,0], [55,0,0,5], [0,40,30,40]],
    #     [[250,400,350], [400,600,350], [200,400,250]],
    #     [[90,75,75,80], [35,85,55,65], [125,95,90,105], [45,110,95,115]],
    #     [[7,2,3],[4,5,6]], #Non square matrix
    #     [[10,19,8,15],[10,18,7,17],[13,16,9,14],[12,19,8,18],[14,17,10,19]],
    #     [[20,22,14,24],[20,19,12,20],[13,10,18,16],[22,23,9,28]],
    #     ]

    # map(termination_test, xs)

    # info("Start mass testing")
    # dimension = [4, 5]
    # sample_size = 10000

    # for testing_dimension in dimension:
    #     for sample in xrange(sample_size):
    #         info("Testing dimension %s. Sample %s" % (testing_dimension, sample))
    #         m = [[random.randrange(100) for i in xrange(testing_dimension)] for j in xrange(testing_dimension)]
    #         termination_test(m)
