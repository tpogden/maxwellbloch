# -*- coding: utf-8 -*-

import nose

def run():
    """
    Run all tests with nose.
    """

    # runs tests in maxwellbloch.tests module
    nose.run(defaultTest="maxwellbloch.tests", argv=['nosetests', '-v'])
