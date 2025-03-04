def initialize():
    '''Initializes the global variables needed for the simulation.
    Note: this function is incomplete, and you may want to modify it'''

    global cur_hedons, cur_health
    global cur_time
    global last_activity, last_activity_duration
    global last_finished
    global bored_with_stars
    global cur_star
    global cur_star_activity
    global last_star_times

    cur_hedons = 0 #User starts at 0 hedons and health points
    cur_health = 0

    cur_star = None
    cur_star_activity = None
    last_star_times = []

    bored_with_stars = False

    last_activity = None
    last_activity_duration = 0

    cur_time = 0

    last_finished = -1000

# Finds if there is an available star that can be used
def star_can_be_taken(activity):
    if cur_star_activity == activity and not bored_with_stars:
        return True
    else:
        return False

# Performs an activity for a duration amount of time
def perform_activity(activity, duration):
    global cur_hedons, cur_health
    global cur_time
    global last_finished, last_activity, last_activity_duration
    global bored_with_stars
    global cur_star_activity
    is_tired = None
    prev = 0

    if duration > 0 and activity in ("running", "textbooks", "resting"):
        if cur_time - last_finished < 120: #check if tired
            is_tired = True
        else:
            is_tired = False

        cur_time += duration #add time

        if activity == "running":              #Code for the activity running
            if last_activity == "running":
                prev = last_activity_duration
                last_activity_duration = duration + prev
            else:
                last_activity_duration = duration

            if duration <= 180 - prev:           #Add Health Points
                cur_health += 3 * duration
            elif prev >= 180:
                cur_health += duration
            elif prev < 180:
                cur_health += (180 - prev) * 3 + (duration - (180 - prev))

            if cur_star_activity != "running" or bored_with_stars:
                if is_tired: #Add/Sub hedons
                    cur_hedons += -2 * duration     #Tired case
                elif duration <= 10:
                    cur_hedons += 2 * duration
                else:
                    cur_hedons += 20 - (duration - 10) * 2
            else: #With Star
                if is_tired: #Add/Sub hedons with Star
                    if duration <= 10:
                        cur_hedons += duration
                    else:
                        cur_hedons += 10 - (duration - 10) * 2
                elif duration <= 10:
                    cur_hedons += 5 * duration
                else:
                    cur_hedons += 50 - (duration - 10) * 2

            cur_star_activity = None #reset current star activity
            last_finished = cur_time
            last_activity = "running"

        elif activity == "textbooks":          #Code for the activity textbooks

            cur_health += 2 * duration         #Add Health Points

            if cur_star_activity != "textbooks" or bored_with_stars:
                if is_tired: #Add/Sub hedons
                    cur_hedons += -2 * duration     #Tired case
                elif duration <= 20:
                    cur_hedons += duration
                else:
                    cur_hedons += 20 - (duration - 20)
            else: #With Star
                if is_tired: #Add/Sub hedons
                    if duration <= 10:
                        cur_hedons += duration
                    else:
                        cur_hedons += 10 - (duration - 10) * 2
                elif duration <= 10:
                    cur_hedons += 4 * duration
                elif duration <= 20:
                    cur_hedons += 40 + (duration - 10)
                else:
                    cur_hedons += 50 - (duration - 20)

            cur_star_activity = None           #reset current star activity
            last_finished = cur_time
            last_activity_duration, last_activity = duration, "textbooks"

        elif activity == "resting":            #Code for the activity resting
            cur_star_activity = None
            last_activity_duration, last_activity = duration, "resting"


    else:
        return None



# Return what current hedons is valued at
def get_cur_hedons():
    global cur_hedons
    return cur_hedons

# Return what current health is valued at
def get_cur_health():
    global cur_health
    return cur_health

# Activates a star for activity
def offer_star(activity):
    global cur_star_activity
    global bored_with_stars
    global last_star_times
    global cur_time
    if activity in ("running", "textbooks", "resting"):
        cur_star_activity = activity
        last_star_times.append(cur_time)
        if len(last_star_times) == 3:
            if last_star_times[2] - last_star_times[0] < 120:
                bored_with_stars = True
                del last_star_times[0]
            else:
                del last_star_times[0]



# Finds the best activity to do at the moment
def most_fun_activity_minute():
    global last_finished, cur_star_activity
    global bored_with_stars
    activities = ["running", "textbooks", "resting"]
    best_activity = None
    hedons = 0  # Calculate hedons for one minute of this activity
    is_tired = None

    if cur_time - last_finished < 120: #check if tired
        is_tired = True
    else:
        is_tired = False

    for activity in activities:

        if activity == "running":
            if cur_star_activity == "running" and not bored_with_stars:
                if is_tired:
                    if hedons < 1:
                        hedons = 1
                        best_activity = "running"
                else:
                    return "running"
            else:
                if is_tired:  # Tired case
                    hedons = -2
                    best_activity = "running"
                else:
                    hedons = 2  # 2 hedons/min otherwise for the first 10 mins
                    best_activity = "running"

        elif activity == "textbooks":
            if cur_star_activity == "textbooks" and not bored_with_stars:
                if is_tired:
                    if hedons < 1:
                        hedons = 1
                        best_activity = "textbooks"
                    elif hedons == 1:
                        best_activity = "running or textbooks"
                else:
                    return "textbooks"
            else:
                if is_tired:
                    if hedons < -2:
                        hedons = -2
                        best_activity = "textbooks"
                else:
                    if hedons < 1:
                        best_activity = "textbooks"
                        hedons = 1

        elif activity == "resting":
            if hedons < 0:
                return "resting"


    return best_activity

#if __name__ == "__main__":
