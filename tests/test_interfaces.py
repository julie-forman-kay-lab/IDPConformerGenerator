"""Test interfaces."""
from datetime import datetime
import pytest

from idpconfgen.core import interfaces as ITF


class TestPrototype:
    """
    Test prototype implementation.

    Following prototype implementation described in:
    https://refactoring.guru/design-patterns/prototype/python/example
    """
    p1 = ITF.Prototype()
    p1.var1 = 1
    p1.var2 = [2, 3, ['4', '6']]
    p1.var3 = datetime.now()
    p1._var4 = 'string'
    p1.circular_reference = ITF.ComponentWithBackReference(p1)
    
    p2 = p1.clone()

    class PTT(ITF.Prototype):
        def __init__(self, varX):
            self.varX = varX
    
    p3 = PTT(datetime.now())
    p3.circular_reference = ITF.ComponentWithBackReference(p3)
    p4 = p3.clone()
    
    @pytest.mark.parametrize(
        'in_p1,in_p2',
        [
            (p1, p2),
            (p1.var2, p2.var2),
            (p1.var3, p2.var3),
            (p1.circular_reference, p2.circular_reference),
            (p1.circular_reference.prototype, p2.circular_reference.prototype),
            (p3, p4),
            (p3.varX, p4.varX),
            ],
        )
    def test_assert_is_not(self, in_p1, in_p2):
        assert in_p1 is not in_p2
    
    @pytest.mark.parametrize(
        'in1,in2',
        [
            (p1, p1),
            (p1._var4, p2._var4),
            ],
        )
    def test_assert_is(self, in1, in2):
        assert in1 is in2
    
    def test_attribute_changes(self):

        self.p1.var1 = 2
        assert not(self.p1.var1 == self.p2.var1)
        
        self.p1.var2[0] = 3
        assert not(self.p1.var2 == self.p2.var2)
        assert self.p1.var2[0] == 3
        assert self.p2.var2[0] == 2

        self.p1.var2[2][0] = 'a'
        assert self.p1.var2[2][0] == 'a'
        assert self.p2.var2[2][0] == '4'
