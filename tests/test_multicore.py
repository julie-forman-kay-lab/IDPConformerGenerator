"""Test libs for multicore operations."""
import os
from functools import partial

import pytest

from idpconfgen.libs import libmulticore as LM


def dummy_generator(a, b=None):
    """Do dummy generator."""
    for i in range(10):
        yield i, a, b


def test_consume_in_list():
    """Test consume generator in list."""
    results = LM.consume_iterable_in_list(dummy_generator, 'foo', b='bar')
    for i, (j, a, b) in enumerate(results):
        assert i == j
        assert a == 'foo'
        assert b == 'bar'


def test_pool_function():
    """Test pool function."""
    execute = LM.pool_function(str, range(10_000), ncores=os.cpu_count())
    result = list(execute)
    assert all(isinstance(i, str) for i in result)
    assert len(result) == 10_000


@pytest.mark.skip(reason="something is halting this")
def test_pool_in_chunks():
    """Test pool in fragments."""
    execute = LM.pool_function_in_chunks(
        str,
        range(10_000),
        ncores=os.cpu_count(),
        chunks=200,
        )
    for chunk in execute:
        assert len(chunk) == 200
        assert isinstance(chunk, list)
        assert isinstance(chunk[0], str)


@pytest.mark.skip(reason="something is halting this")
def test_pool_in_chunks_nested():
    """Test pool in fragment from nested result."""
    execute = partial(LM.consume_iterable_in_list, dummy_generator)
    execute_pool = LM.pool_function_in_chunks(
        execute,
        range(20),
        ncores=os.cpu_count(),
        chunks=5,
        )

    for chunk in execute_pool:

        # fragment of results
        assert len(chunk) == 5
        assert isinstance(chunk, list)

        # this is a list containing the results from `dummy_generator`
        # list of tuples, contains 10 elements because `dummy_generator`
        # yields 10x
        assert isinstance(chunk[0], list)
        assert len(chunk[0]) == 10

        # each yielded results is a tuple of 3 elements
        assert isinstance(chunk[0][0], tuple)
        assert len(chunk[0][0]) == 3


@pytest.mark.skip(reason="something is halting this")
def test_pool_in_chunks_flatten():
    """Test pool in fragment from nested result."""
    execute = partial(LM.consume_iterable_in_list, dummy_generator)
    execute_pool = partial(
        LM.pool_function_in_chunks,
        execute,
        range(20),
        ncores=os.cpu_count(),
        chunks=5,
        )

    def assrt(flat):
        for result in flat:
            assert isinstance(result, tuple)

    LM.flat_results_from_chunk(execute_pool, assrt)
