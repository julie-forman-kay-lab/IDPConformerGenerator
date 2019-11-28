"""Objects interfaces."""


class _ComponentWithBackReference:
    """
    https://refactoring.guru/design-patterns/prototype/python/example
    """
    def __init__(self, protopype):
        self._prototype = prototype

    @property
    def prototype(self):
        return self._prototype

    @prototype.setter
    def prototype(self, value):
        self._prototype = value


class Prototype:
    """Prototype interface."""
    def __init__(self):
        self._circular_reference = None

    @property
    def circular_reference(self):
        return self._circular_reference

    @property
    def circular_reference(self, value):
        assert isinstance(value, _ComponentWithBackReference)
        self._circular_reference = value
