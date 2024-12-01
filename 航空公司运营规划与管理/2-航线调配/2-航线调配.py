import datetime
import numpy as np

FlightMatrix = [
    [125, 'JFK', '07:25', 'SFO', '09:55', 5.5],
    [110, 'ATL', '08:10', 'JFK', '10:40', 2.5],
    [113, 'MIA', '09:10', 'JFK', '12:10', 3],
    [131, 'JFK', '09:30', 'ATL', '12:00', 2.5],
    [105, 'SFO', '09:50', 'JFK', '18:20', 5.5],
    [138, 'JFK', '12:30', 'BOS', '14:00', 1.5],
    [111, 'ATL', '13:10', 'JFK', '15:40', 2.5],
    [114, 'MIA', '14:30', 'JFK', '17:30', 3],
    [118, 'BOS', '15:00', 'JFK', '16:30', 1.5],
    [135, 'JFK', '15:10', 'MIA', '18:10', 3],
    [133, 'JFK', '18:05', 'ATL', '20:35', 2.5],
    [136, 'JFK', '18:10', 'MIA', '21:10', 3]
]

DepartureMinutes = np.zeros(len(FlightMatrix))
ArrivalMinutes = np.zeros(len(FlightMatrix))

for i in range(len(FlightMatrix)):
    dep_time = datetime.datetime.strptime(FlightMatrix[i][2], '%H:%M')
    DepartureMinutes[i] = dep_time.hour * 60 + dep_time.minute

    arr_time = datetime.datetime.strptime(FlightMatrix[i][4], '%H:%M')
    ArrivalMinutes[i] = arr_time.hour * 60 + arr_time.minute


