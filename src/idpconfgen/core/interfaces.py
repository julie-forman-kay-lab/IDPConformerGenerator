"""Objects interfaces."""
from copy import deepcopy


class ComponentWithBackReference:
    """
    https://refactoring.guru/design-patterns/prototype/python/example
    """
    def __init__(self, prototype):
        self._prototype = prototype

    @property
    def prototype(self):
        return self._prototype

    @prototype.setter
    def prototype(self, value):
        self._prototype = value


class Prototype:
    """
    Prototype interface.

    See example:
        https://refactoring.guru/design-patterns/prototype/python/example
    """
    def __init__(self):
        self._circular_reference = None

    @property
    def circular_reference(self):
        return self._circular_reference

    @circular_reference.setter
    def circular_reference(self, value):
        assert isinstance(value, ComponentWithBackReference)
        self._circular_reference = value

    def clone(self):
        """
        Clones object.
        """
        self.circular_reference = deepcopy(self.circular_reference)
        self.circular_reference.prototype = self
        return deepcopy(self)
