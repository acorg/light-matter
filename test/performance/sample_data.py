# Sample performance test results.

RESULT1 = {
    "UTC": "2015-01-09 16:02:07",
    "startTestRunTime": 1420819327.9003059864,
    "description": "<not given>",
    "results": {
        "performance.perf_database.TestDatabase.testAddSubjects": {
            "status": "success",
            "details": 0.0606100559,
            "elapsed": 0.0608661175,
        },
        "performance.perf_database.TestDatabase.testChecksumEmpty": {
            "status": "success",
            "details": 0.0000140667,
            "elapsed": 0.000054121,
        },
        "performance.perf_findSelf.TestFindSelf.testFindIdentical": {
            "status": "success",
            "details": {
                "score": 48,
                "result": True,
            },
            "elapsed": 0.001333952
        },
        "performance.perf_database.TestDatabase.testCreation": {
            "status": "success",
            "elapsed": 0.000043869
        },
        "performance.perf_database.TestDatabase.testChecksum": {
            "status": "success",
            "details": 0.0157001019,
            "elapsed": 0.0762839317,
        },
    },
    "elapsed": 0.1387579441,
    "testCount": 5
}

RESULT2 = {
    "UTC": "2015-01-09 17:02:07",
    "startTestRunTime": 1420819327.9003059864,
    "description": "<not given>",
    "results": {
        "performance.perf_database.TestDatabase.testAddSubjects": {
            "status": "success",
            "details": 0.1606100559,
            "elapsed": 0.1608661175,
        },
        "performance.perf_database.TestDatabase.testChecksumEmpty": {
            "status": "success",
            "details": 0.1000140667,
            "elapsed": 0.100054121,
        },
        "performance.perf_findSelf.TestFindSelf.testFindIdentical": {
            "status": "success",
            "details": {
                "score": 48,
                "result": True,
            },
            "elapsed": 0.101333952
        },
        "performance.perf_database.TestDatabase.testCreation": {
            "status": "success",
            "elapsed": 0.100043869
        },
        "performance.perf_database.TestDatabase.testChecksum": {
            "status": "success",
            "details": 0.1157001019,
            "elapsed": 0.1762839317,
        },
        "performance.perf_database.TestDatabase.testSomething": {
            "status": "success",
            "details": 0.1157001019,
            "elapsed": 0.1762839317,
        },
    },
    "elapsed": 0.2387579441,
    "testCount": 5
}
