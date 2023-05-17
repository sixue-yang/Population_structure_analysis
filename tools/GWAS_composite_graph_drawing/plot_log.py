# _author_ = 'sixueyang'
# _date_ = 2023/5/16 14:21
import datetime



class Log:
    def __init__(self, level='INFO', filename='log.txt', format='{time} {level}: {message}'):
        self.level = level.upper()
        self.format = format
        self.closed = False

        # open log file
        self.log_file = open(filename, 'a')

    def log(self, message, level=None):
        if self.closed:
            raise ValueError("Log file is closed")

        time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        level = level.upper() if level else self.level

        log_message = self.format.format(time=time, level=level, message=message)

        if level == 'DEBUG':
            print(log_message)
        elif level == 'INFO' or level == 'WARNING':
            print(log_message)
        elif level == 'ERROR' or level == 'CRITICAL':
            print(log_message)

        # write log message to file
        print(log_message,file=self.log_file)
        self.log_file.flush()
        # self.log_file.write(log_message + '\n')

    def debug(self, message):
        self.log(message, 'DEBUG')

    def info(self, message):
        self.log(message, 'INFO')

    def warning(self, message):
        self.log(message, 'WARNING')

    def error(self, message):
        self.log(message, 'ERROR')

    def critical(self, message):
        self.log(message, 'CRITICAL')

    def close(self):
        self.closed = True
        # close log file
        self.log_file.close()



