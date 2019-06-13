"""
    pyexcel.io
    ~~~~~~~~~~~~~~~~~~~

    Uniform interface for import/export data from different db interfaces

    :copyright: (c) 2014-2015 by Onni Software Ltd.
    :license: New BSD License, see LICENSE for more details
"""

class DjangoModel:
    pass

class AlchemyModel:
    pass

EXPORT_HANDLERS = {
    "django": DjangoModel,
    "sqlalchemy": AlchemyModel
}