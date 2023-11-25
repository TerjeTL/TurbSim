from PyQt6.QtWidgets import *

# Only needed for access to command line arguments
import sys

# You need one (and only one) QApplication instance per application.
# Pass in sys.argv to allow command line arguments for your app.
# If you know you won't use command line arguments QApplication([]) works too.
app = QApplication(sys.argv)

# Subclass QMainWindow to customize your application's main window
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Bezier-C")
        button = QPushButton("Press Me!")
        button.setCheckable = True
        button.clicked.connect(self.on_button_click)

        # Set the central widget of the Window.
        self.setCentralWidget(button)
    
    def on_button_click(self):
        print("Clicked")

# Create a Qt widget, which will be our window.
window = MainWindow()
window.show()  # IMPORTANT!!!!! Windows are hidden by default.

# Start the event loop.
app.exec()


# Your application won't reach here until you exit and the event
# loop has stopped.
